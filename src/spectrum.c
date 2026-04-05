//#define SPECTRUM_DEBUG 1

// Spectral analysis service - far from complete - for ka9q-radio's radiod
// Copyright 2023-2025, Phil Karn, KA9Q
#define _GNU_SOURCE 1
#include <assert.h>
#include <pthread.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "misc.h"
#include "iir.h"
#include "filter.h"
#include "radio.h"
#include "window.h"

//#define RICE 1
//#define SPECTRUM_CLIP 1
//#define FIXED_STEP  1

static void generate_window(struct channel *);
static void setup_real_fft(struct channel *);
static void setup_complex_fft(struct channel *);
static void setup_wideband(struct channel *);
static void setup_narrowband(struct channel *);
static void narrowband_poll(struct channel *);
static void wideband_poll(struct channel *);
static int wideband_poll_n(struct channel *, int count, double gain);
static void wideband_update_scaling(struct channel *);
#if RICE
static void rice(struct channel *);
#endif

// Spectrum analysis thread
int demod_spectrum(void *arg){
  struct channel * const chan = arg;
  assert(chan != NULL);
  if(chan == NULL)
    return -1;

  pthread_mutex_init(&chan->status.lock,NULL);
  pthread_mutex_lock(&chan->status.lock);
  chan->status.output_interval = 0; // No automatic status updates
  chan->status.output_timer = 0; // No automatic status updates
  chan->output.silent = true; // we don't send anything there
  pthread_mutex_unlock(&chan->status.lock);

  // Parameters set by system input side
  assert(Blocktime != 0);
  realtime(chan->prio - 10); // Drop below demods
  if(chan->spectrum.fft_avg <= 0)
    chan->spectrum.fft_avg = 1;     // force legal

  bool restart_needed = false;
  bool response_needed = true;
  // Watch for parameter changes and do them in the loop so we don't have to force a restart
  enum window_type window_type = -1; // force generation on first loop
  double rbw = -1;
  int bin_count = -1;
  int crossover = -1;
  double shape = -1;
  int timeout = 0;

  // Main loop
  do {
    response(chan,response_needed);
    response_needed = false;

    // Look on the single-entry command queue, grab it atomically and execute it
    pthread_mutex_lock(&chan->status.lock);
    // Look on the command queue and grab just one atomically
    for(int i=0;i < CQLEN; i++){
      if(chan->commands[i].buffer != NULL){
	restart_needed = decode_radio_commands(chan,chan->commands[i].buffer,
					       chan->commands[i].length);
	FREE(chan->commands[i].buffer);
	chan->commands[i].length = 0;
	response_needed = true;
	break;
      }
    }
    pthread_mutex_unlock(&chan->status.lock);

    // Must handle possible parameter changes from decode_radio_commands() BEFORE executing the downconverter,
    // which will act immediately on those changes. Otherwise segfaults occur when crossing between wideband and narrowband
    // modes because things are not properly set up for the poll when it comes
    if((chan->spectrum.rbw > chan->spectrum.crossover) != (rbw > crossover)) // note nested booleans
      chan->spectrum.fft_n = -1; // force setup

    if(chan->spectrum.rbw != rbw || chan->spectrum.bin_count != bin_count)
      chan->spectrum.fft_n = -1; // force setup;

    // fairly major reinitialization required
    if(chan->spectrum.fft_n <= 0){
      chan->spectrum.accum_remaining = 0; // cancel any in-progress accumulation
      if(chan->spectrum.plan != NULL)
	fftwf_destroy_plan(chan->spectrum.plan); // will be regenerated on first poll
      chan->spectrum.plan = NULL;
      FREE(chan->spectrum.window); // force regeneration on first poll
      if(chan->spectrum.rbw > chan->spectrum.crossover)
	setup_wideband(chan);
      else
	setup_narrowband(chan);
    } else if(chan->spectrum.window_type != window_type
	      || (chan->spectrum.shape != shape && (chan->spectrum.window_type == KAISER_WINDOW
						    || chan->spectrum.window_type == GAUSSIAN_WINDOW))){
      FREE(chan->spectrum.window); // force regeneration
    }
    // End of parameter checking and (re)initialization
    if(restart_needed)
      break; // restart or terminate
    int r = downconvert(chan);
    if(r == -1)
      break; // restart needed
    if(r == 1 && !response_needed && chan->spectrum.accum_remaining <= 0)
      continue; // channel inactive, no pending response, no accumulation in progress

    // r == 0 is normal return
    // Process receiver data only in narrowband mode
    if(r == 0 && chan->spectrum.rbw <= chan->spectrum.crossover && chan->baseband != NULL){
      if(chan->spectrum.ring == NULL || chan->spectrum.ring_size < chan->spectrum.fft_avg * chan->spectrum.fft_n){

	// Need a new or bigger baseband ring buffer
	if(chan->spectrum.ring == NULL)
	  chan->spectrum.ring_idx = 0; // ? no need to reset on growth vs start?

	chan->spectrum.ring_size = chan->spectrum.fft_avg * chan->spectrum.fft_n;
	assert(chan->spectrum.ring_size > 0);
	void *old = chan->spectrum.ring;
	int old_ring_size = chan->spectrum.ring_size;
	chan->spectrum.ring = realloc(chan->spectrum.ring, chan->spectrum.ring_size * sizeof *chan->spectrum.ring);
	if(chan->spectrum.ring == NULL)
	  FREE(old); // emulate reallocf() in case the realloc fails, though we'll crash anyway
	// Clear the new space to avoid display glitches
	memset(chan->spectrum.ring + old_ring_size, 0, (chan->spectrum.ring_size - old_ring_size) * sizeof *chan->spectrum.ring);
      }
      assert(chan->spectrum.ring != NULL);
      for(int i=0; i < chan->sampcount; i++){
	chan->spectrum.ring[chan->spectrum.ring_idx++] = chan->baseband[i];
	if(chan->spectrum.ring_idx == chan->spectrum.ring_size)
	  chan->spectrum.ring_idx = 0; // wrap around
      }
      timeout -= chan->sampcount;
      if(timeout < 0)
	timeout = 0;
    }
    // Ensure bin_data buffer exists and is correctly sized
    if(chan->spectrum.bin_data == NULL || chan->spectrum.bin_count != bin_count){
      void *old = chan->spectrum.bin_data;
      chan->spectrum.bin_data = realloc(chan->spectrum.bin_data, chan->spectrum.bin_count * sizeof *chan->spectrum.bin_data);
      if(chan->spectrum.bin_data == NULL)
	FREE(old); // emulate reallocf()
    }

    if(response_needed && chan->spectrum.rbw > chan->spectrum.crossover){
      // Wideband mode: start incremental accumulation across frames
      int const fft_avg = chan->spectrum.fft_avg <= 0 ? 1 : chan->spectrum.fft_avg;
      memset(chan->spectrum.bin_data, 0, chan->spectrum.bin_count * sizeof *chan->spectrum.bin_data);
      chan->spectrum.accum_remaining = fft_avg;
      // Precompute loop constants once here; fft_n/overlap/samprate don't change mid-accumulation
      // (a parameter change forces fft_n=-1 which cancels accum_remaining first)
      // fft_step = samples advanced per FFT frame = fft_n * (1 - overlap)
      // e.g. with 50% overlap and fft_n=1250: step=625, so 2 frames per 1250-sample block
      chan->spectrum.accum_step = lrint(chan->spectrum.fft_n * (1. - chan->spectrum.overlap));
      int const block_samples = lrint(chan->frontend->samprate * Blocktime);
      int const batch_max = block_samples / chan->spectrum.accum_step;
      chan->spectrum.accum_batch_max = batch_max < 1 ? 1 : batch_max;
      // gain normalizes the accumulated sum: divide by total frames and by fft_n^2
      // (fft_n^2 because cnrm gives |X|^2 which scales as N^2 for a normalized window)
      chan->spectrum.accum_gain = 1./(double)((int64_t)fft_avg * chan->spectrum.fft_n * chan->spectrum.fft_n);
      response_needed = false; // defer until accumulation completes
    }

    if(chan->spectrum.accum_remaining > 0 && chan->spectrum.rbw > chan->spectrum.crossover){
      // Wideband incremental: process a batch of FFT frames from the latest data.
      // We spread the total fft_avg frames across multiple downconvert cycles so that
      // each batch reads genuinely new data from the A/D ring buffer rather than
      // re-reading the same samples from a stale snapshot.
      int batch = chan->spectrum.accum_batch_max;
      if(batch > chan->spectrum.accum_remaining)
	batch = chan->spectrum.accum_remaining;
      chan->spectrum.accum_remaining -= wideband_poll_n(chan, batch, chan->spectrum.accum_gain);
      if(chan->spectrum.accum_remaining <= 0){
	wideband_update_scaling(chan);
	response_needed = true; // accumulation complete, send response
#ifdef RICE
	rice(chan);
#endif
      }
    } else if(response_needed){
      // Narrowband mode or immediate response needed
      if(chan->spectrum.rbw <= chan->spectrum.crossover){
	narrowband_poll(chan);
      } else {
	wideband_poll(chan); // fallback: single-shot
      }
#ifdef RICE
      rice(chan);
#endif
    }
    // Remember new values in case they change next time
    rbw = chan->spectrum.rbw;
    bin_count = chan->spectrum.bin_count;
    crossover = chan->spectrum.crossover;
    window_type = chan->spectrum.window_type;
    shape = chan->spectrum.shape;
  } while(true);

  if(Verbose > 1)
    fprintf(stderr,"%s exiting\n",chan->name);

  chan->spectrum.fft_n = 0;
  delete_filter_output(&chan->filter.out);
  chan->baseband = NULL;
  if(chan->spectrum.plan)
    fftwf_destroy_plan(chan->spectrum.plan);
  chan->spectrum.plan = NULL;
  FREE(chan->spectrum.window);
  fftwf_free(chan->spectrum.fft_in_r);
  chan->spectrum.fft_in_r = NULL;
  fftwf_free(chan->spectrum.fft_in_c);
  chan->spectrum.fft_in_c = NULL;
  fftwf_free(chan->spectrum.fft_out);
  chan->spectrum.fft_out = NULL;
  for(int i=0; i < CQLEN; i++){
    FREE(chan->commands[i].buffer);
    chan->commands[i].length = 0;
  }
  FREE(chan->spectrum.bin_data);
  FREE(chan->spectrum.ring);
  chan->spectrum.ring_size = 0;
  return 0;
}

static void narrowband_poll(struct channel *chan){
  // Narrowband mode poll

  if(chan->spectrum.ring == NULL)
      return; // Needed

  struct frontend const * restrict const frontend = chan->frontend;
  if(frontend == NULL)
    return;

  int const bin_count = chan->spectrum.bin_count;
  float * restrict const bin_data = chan->spectrum.bin_data;
  assert(bin_data != NULL); // allocated just before we're called
  memset(bin_data,0, bin_count * sizeof *bin_data); // zero output data

  if(chan->spectrum.plan == NULL)
    setup_complex_fft(chan); // narrowband always uses complex

  fftwf_plan restrict const plan = chan->spectrum.plan;
  assert(plan != NULL);

  if(chan->spectrum.window == NULL)
    generate_window(chan);

  float const * restrict const window = chan->spectrum.window;

  // Most recent data from receive ring buffer
  float complex const * restrict const ring = chan->spectrum.ring;
  int const ring_size = chan->spectrum.ring_size;

  int const fft_n = chan->spectrum.fft_n;
  assert(fft_n > 0); // should be set by narrowband_setup()

  int rp = chan->spectrum.ring_idx - fft_n;
  if(rp < 0)
    rp += ring_size;

  // Use persistent buffers allocated in setup_narrowband; avoids fftwf_alloc/free every call
  float complex * restrict fft_in = chan->spectrum.fft_in_c;
  float complex * restrict fft_out = chan->spectrum.fft_out;
  assert(fft_in != NULL && fft_out != NULL);

  int fft_avg = chan->spectrum.fft_avg;
  fft_avg = fft_avg <= 0 ? 1 : fft_avg; // force it valid

  // scale each bin value for our FFT
  // squared because the we're scaling the output of complex norm, not the input bin values
  // Unlike wideband, no adjustment for a real front end because the downconverter corrects the gain
  double const gain = 1.0 / ((double)fft_n * fft_n * fft_avg);
  // rp_step: how far back to move rp between iterations to achieve the configured overlap
  int const rp_step = lrint(fft_n * (2. - chan->spectrum.overlap)); // hoisted out of loop

  for(int iter=0; iter < fft_avg; iter++){
    // Copy and window raw baseband
    for(int i = 0; i < fft_n; i++){
      assert(rp >= 0 && rp < ring_size);
      fft_in[i] = ring[rp++] * window[i];
      if(rp >= ring_size)
	rp -= ring_size;
    }
    fftwf_execute_dft(plan,fft_in,fft_out);
    // DC to Nyquist-1, then -Nyquist to -1
    int fr = 0;
    for(int i=0; i < bin_count; i++){
      if(i == bin_count/2)
	fr = fft_n - i; // skip over excess FFT bins at edges
      assert(fr >= 0 && fr < fft_n);
      double const p = cnrm((double complex)fft_out[fr++]); // use double for improved accuracy when summing?
      if(isfinite(p))
	bin_data[i] += gain * p; // Don't pollute with infinities or NANs
    }
    // rp now points to *next* buffer, so move it back between 1 and 2 buffers depending on overlap
    rp -= rp_step;
    if(rp < 0)
      rp += ring_size;
  }
  double min_power = INFINITY;
  double max_power = 0;

  // scaling for byte bin format
  for(int i=0; i < bin_count; i++){
    if(bin_data[i] < min_power)
      min_power = bin_data[i];
    if(bin_data[i] > max_power)
      max_power = bin_data[i];
  }
  if(max_power > 0 && min_power > 0){
#if SPECTRUM_CLIP
    chan->spectrum.base = power2dB(chan->sig.n0 * chan->spectrum.noise_bw);
#else
    chan->spectrum.base = power2dB(min_power);
#endif
#if FIXED_STEP
    chan->spectrum.step = 0.5; // 0.5 dB fixed
#else
    chan->spectrum.step = (1./256.) * (power2dB(max_power) - chan->spectrum.base); // dB range
#endif
  }
}

// Process 'count' FFT frames from the most recent data in the frontend A/D buffer
// Reads backward from write pointer with overlap between frames
// Adds results to bin_data (caller must zero bin_data before first call of a new accumulation)
// gain should include the 1/total_avg normalization factor
// Returns number of frames actually processed
static int wideband_poll_n(struct channel *chan, int count, double gain){
  if(chan == NULL || count <= 0)
    return 0;

  struct frontend const * restrict const frontend = chan->frontend;
  if(frontend == NULL)
    return 0;

  // These can happen if we're called too early, before allocations
  struct filter_in const * const restrict master = chan->filter.out.master;
  if(master == NULL)
    return 0;

  int const fft_n = chan->spectrum.fft_n;
  assert(fft_n > 0); // should be set by setup_wideband()

  if(chan->spectrum.plan == NULL){
    if(chan->frontend->isreal)
      setup_real_fft(chan);
    else
      setup_complex_fft(chan);
  }
  fftwf_plan restrict const plan = chan->spectrum.plan;
  assert(plan != NULL);

  // should be set up just before we were called
  int const bin_count = chan->spectrum.bin_count;
  float * restrict const bin_data = chan->spectrum.bin_data;
  assert(bin_count > 0 && bin_data != NULL);

  if(chan->spectrum.window == NULL)
    generate_window(chan);
  float const * restrict const window = chan->spectrum.window;
  assert(window != NULL);

  // Asynchronously read newest data from input buffer
  // Look back from the most recent write pointer to allow room for overlapping windows
  // scale fft bin shift down to size of analysis FFT, which is smaller than the input FFT
  int shift = (int)(chan->filter.bin_shift * (int64_t)fft_n / master->points);
  int processed = 0;

  if(frontend->isreal){
    // Point into raw SDR A/D input ring buffer
    // We're reading from a mirrored buffer so it will automatically wrap back to the beginning
    // as long as it doesn't go past twice the buffer length
    float const *input = frontend->in.input_write_pointer.r - fft_n; // 1 FFT buffer back
    if(input < (float *)frontend->in.input_buffer)
      input += frontend->in.input_buffer_size / sizeof *input; // wrap backward

    // Use persistent buffers allocated in setup_wideband; avoids fftwf_alloc/free every call
    float * restrict fft_in = chan->spectrum.fft_in_r;
    float complex * restrict fft_out = chan->spectrum.fft_out;
    assert(fft_in != NULL && fft_out != NULL);
    // scale each bin value for our FFT
    // squared because we're scaling the output of complex norm, not the input bin values
    // +3dB to include the virtual conjugate spectrum (real FFT discards the negative frequencies)
    double const real_gain = 2.0 * gain;
    int const fft_step = lrint(fft_n * (1. - chan->spectrum.overlap)); // samples per step, hoisted out of loop

    for(int iter = 0; iter < count; iter++){
      // Copy and window raw A/D
      if(shift >= 0){
	// Upright spectrum
	for(int i=0; i < fft_n; i++)
	  fft_in[i] = window[i] * input[i];
      } else {
	// Invert spectrum by flipping sign of every other sample - UNTESTED
	// equivalent to multiplication by a sinusoid at the Nyquist rate
	// If FFT N is odd, just forget the odd last sample.
	// We don't have to track the sign flip phase because we're only summing energy
	for(int i=0; i < fft_n-1; i += 2){
	  fft_in[i] = window[i] * input[i];
	  fft_in[i+1] = -window[i+1] * input[i+1];
	}
	if(fft_n & 1)
	  fft_in[fft_n-1] = 0;
	shift = -shift; // make positive
      }
      fftwf_execute_dft_r2c(plan,fft_in,fft_out);

      // Spectrum is always right side up so shift is never negative
      // Start with DC + positive frequencies, then wrap to negative
      int binp = shift;
      assert(binp >= 0);
      for(int i=0;i < bin_count && binp < fft_n/2+1 ; i++,binp++){
	if(i == bin_count/2)
	  binp -= bin_count; // crossed into negative output range, wrap input back to lowest frequency requested
	double const p = cnrm(fft_out[binp]);
	if(isfinite(p))
	  bin_data[i] += real_gain * p;
      }
      processed++;
      input -= fft_step; // move back by one step (accounts for overlap)
      if(input < (float *)frontend->in.input_buffer)
	input += frontend->in.input_buffer_size / sizeof *input; // wrap backward
    }
  } else {
    // Complex front end (frontend->isreal == false) - UNTESTED
    // Find starting point to read in input A/D stream
    float complex const * restrict input = frontend->in.input_write_pointer.c - fft_n; // 1 buffer back
    input += (input < (float complex *)frontend->in.input_buffer) ? frontend->in.input_buffer_size / sizeof *input : 0; // backward wrap

    // Use persistent buffers allocated in setup_wideband; avoids fftwf_alloc/free every call
    float complex * restrict fft_in = chan->spectrum.fft_in_c;
    float complex * restrict fft_out = chan->spectrum.fft_out;
    assert(fft_in != NULL && fft_out != NULL);
    int const fft_step = lrint(fft_n * (1. - chan->spectrum.overlap)); // samples per step, hoisted out of loop

    for(int iter = 0; iter < count; iter++){
      // Copy and window raw A/D
      for(int i=0; i < fft_n; i++)
	fft_in[i] = window[i] * input[i];

      fftwf_execute_dft(plan,fft_in,fft_out);

      // Copy requested bins to output, starting with requested center frequency
      // FFTW complex output layout: [0]=DC, [1..N/2-1]=positive, [N/2]=Nyquist, [N/2+1..N-1]=negative
      // bin_data layout: [0..bin_count/2-1]=positive freqs, [bin_count/2..bin_count-1]=negative freqs
      // shift maps the channel's bin_shift (in master FFT bins) to this smaller analysis FFT
      int binp;
      int i = 0;

      if(shift >= 0){
	// Starts in positive spectrum
	if(shift >= fft_n/2)
	  goto skip; // starts past end of spectrum, nothing to return
	binp = shift;
      } else {
	// shift < 0, starts in negative spectrum
	if(-shift >= fft_n)
	  goto skip; // nothing overlaps, quit
	if(-shift >= fft_n/2){ // before start of input spectrum
	  i = -shift - fft_n/2;
	  binp = fft_n/2; // start input at lowest negative frequency
	} else {
	  binp = fft_n + shift;
	}
      }
      for(; i < bin_count; i++){
	if(i == bin_count/2){
	  // Crossed midpoint: output switches from positive to negative frequencies.
	  // In FFTW layout negative freqs are at the top of the array (indices N/2+1..N-1),
	  // so jump binp back by bin_count to land at the right negative-frequency bin.
	  binp -= bin_count;
	  if(binp < 0)
	    binp += fft_n;
	}
	assert(binp >= 0 && binp < fft_n && i >= 0 && i < bin_count);
	double const p = cnrm(fft_out[binp]);
	if(isfinite(p))
	  bin_data[i] += gain * p;

	if(++binp == fft_n)
	  binp = 0; // wrap to DC
      }
    skip:;
      processed++;
      // Back to previous buffer
      input -= fft_step; // move back by one step (accounts for overlap)
      if(input < (float complex *)frontend->in.input_buffer)
	input += frontend->in.input_buffer_size / sizeof *input; // backward wrap
    }
  }
  return processed;
}

// Original single-shot wrapper: process all fft_avg frames at once from existing buffer.
// Used as a fallback when there is no pending incremental accumulation.
static void wideband_poll(struct channel *chan){
  if(chan == NULL)
    return;

  // should be set up just before we were called
  int const bin_count = chan->spectrum.bin_count;
  float * restrict const bin_data = chan->spectrum.bin_data;
  if(bin_count <= 0 || bin_data == NULL)
    return;
  memset(bin_data,0, bin_count * sizeof *bin_data); // zero output data

  int const fft_avg = chan->spectrum.fft_avg <= 0 ? 1 : chan->spectrum.fft_avg; // force it valid
  int const fft_n = chan->spectrum.fft_n;
  // scale each bin value for our FFT
  // squared because we're scaling the output of complex norm, not the input bin values
  double const gain = 1./(double)((int64_t)fft_avg * fft_n * fft_n);

  wideband_poll_n(chan, fft_avg, gain);
  wideband_update_scaling(chan);
}

// Update base/step scaling after accumulation is complete
static void wideband_update_scaling(struct channel *chan){
  int const bin_count = chan->spectrum.bin_count;
  float const * restrict const bin_data = chan->spectrum.bin_data;

  double min_power = INFINITY;
  double max_power = 0;

  for(int i=0; i < bin_count; i++){
    if(bin_data[i] < min_power)
      min_power = bin_data[i];
    if(bin_data[i] > max_power)
      max_power = bin_data[i];
  }
  if(max_power > 0 && min_power > 0){
#if SPECTRUM_CLIP
    chan->spectrum.base = power2dB(chan->sig.n0 * chan->spectrum.noise_bw);
#else
    chan->spectrum.base = power2dB(min_power);
#endif
#if FIXED_STEP
    chan->spectrum.step = 0.5;
#else
    chan->spectrum.step = ldexp(power2dB(max_power) - chan->spectrum.base,-8);
#endif
  }
}

// Fill a buffer with compact frequency bin data, 1 byte each
// Each unsigned byte in the output gives the bin power above 'base' decibels, in steps of 'step' dB
// Unlike the regular float32 format, these bins run uniformly from lowest frequency to highest frequency
void encode_byte_data(struct channel const *chan, uint8_t *buffer){
  assert(chan != NULL && buffer != NULL);

  int const bin_count = chan->spectrum.bin_count;
  double scale = 1./chan->spectrum.step;

  int wbin = bin_count/2; // nyquist freq is most negative
  for(int i=0; i < bin_count; i++){
    double x = scale * (power2dB(chan->spectrum.bin_data[wbin++]) - chan->spectrum.base);
    if(x < 0)
      x = 0;
    else if(x > 255)
      x = 255;

    buffer[i] = (uint8_t)llrint(x);
    if(wbin == bin_count)
      wbin = 0;  // Continuing through dc and positive frequencies
  }
}
// Generate normalized sampling window
// the generation functions are symmetric so lengthen them by one point to make them periodic
static void generate_window(struct channel *chan){
  FREE(chan->spectrum.window);

  chan->spectrum.window = malloc((1 + chan->spectrum.fft_n) * sizeof *chan->spectrum.window);
  assert(chan->spectrum.window != NULL);
  switch(chan->spectrum.window_type){
  default:
  case KAISER_WINDOW: // If β == 0, same as rectangular
    make_kaiserf(chan->spectrum.window,chan->spectrum.fft_n+1,chan->spectrum.shape);
    break;
  case RECT_WINDOW:
    for(int i=0; i < chan->spectrum.fft_n; i++)
      chan->spectrum.window[i] = 1;
    break;
  case BLACKMAN_WINDOW:
    for(int i=0; i < chan->spectrum.fft_n; i++)
      chan->spectrum.window[i] = blackman_window(i,chan->spectrum.fft_n+1);
    break;
  case EXACT_BLACKMAN_WINDOW:
    for(int i=0; i < chan->spectrum.fft_n; i++)
      chan->spectrum.window[i] = exact_blackman_window(i,chan->spectrum.fft_n+1);
    break;
  case BLACKMAN_HARRIS_WINDOW:
    for(int i=0; i < chan->spectrum.fft_n; i++)
      chan->spectrum.window[i] = blackman_harris_window(i,chan->spectrum.fft_n+1);
    break;
  case GAUSSIAN_WINDOW:
    // Reuse kaiser β as σ parameter
    // note σ = 0 is a pathological value for gaussian, it's an impulse with infinite sidelobes
    gaussian_window_alpha(chan->spectrum.window, chan->spectrum.fft_n+1,chan->spectrum.shape, false); // we normalize them all below
    break;
  case HANN_WINDOW:
    for(int i=0; i < chan->spectrum.fft_n; i++)
      chan->spectrum.window[i] = hann_window(i,chan->spectrum.fft_n+1);
    break;
  case HAMMING_WINDOW:
    for(int i=0; i < chan->spectrum.fft_n; i++)
      chan->spectrum.window[i] = hamming_window(i,chan->spectrum.fft_n+1);
    break;
  case HP5FT_WINDOW:
    for(int i=0; i < chan->spectrum.fft_n; i++)
      chan->spectrum.window[i] = hp5ft_window(i,chan->spectrum.fft_n+1);
    break;
  }
  normalize_windowf(chan->spectrum.window,chan->spectrum.fft_n);

  // Compute noise bandwidth of each bin in bins
  chan->spectrum.noise_bw = 0;
  for(int i=0; i < chan->spectrum.fft_n; i++)
    chan->spectrum.noise_bw += (double)chan->spectrum.window[i] * chan->spectrum.window[i];

  // Scale to the actual bin bandwidth
  // This also has to be divided by the square of the sum of the window values, but that's already normalized to 1
  chan->spectrum.noise_bw *= chan->spectrum.rbw / chan->spectrum.fft_n;
}

static void setup_wideband(struct channel *chan){
  // Direct Wideband mode. Setup FFT to work on raw A/D input
  // What can we do about unfriendly sizes? Anything?
  chan->spectrum.fft_n = lrint(chan->frontend->samprate / chan->spectrum.rbw);
  chan->output.samprate = 0; // Not meaningful
  chan->output.channels = 0;
  if(Verbose > 1)
    fprintf(stderr,"%s wide spectrum: center %'.3lf Hz bin count %u, rbw %.1lf Hz, samprate %u Hz fft size %u\n",
	    chan->name,chan->tune.freq,chan->spectrum.bin_count,chan->spectrum.rbw,chan->output.samprate,chan->spectrum.fft_n);

  FREE(chan->spectrum.ring); // not needed
  chan->spectrum.ring_size = 0;
  // Dummy just so downconvert() will block on each frame
  delete_filter_output(&chan->filter.out);
  int r = create_filter_output(&chan->filter.out,&chan->frontend->in,NULL,0,SPECTRUM);
  assert(r == 0);
  (void)r;
  // Allocate persistent FFT I/O buffers sized for this fft_n; avoids per-call heap churn
  fftwf_free(chan->spectrum.fft_in_r);
  chan->spectrum.fft_in_r = NULL;
  fftwf_free(chan->spectrum.fft_in_c);
  chan->spectrum.fft_in_c = NULL;
  fftwf_free(chan->spectrum.fft_out);
  if(chan->frontend->isreal)
    chan->spectrum.fft_in_r = fftwf_alloc_real(chan->spectrum.fft_n);
  else
    chan->spectrum.fft_in_c = fftwf_alloc_complex(chan->spectrum.fft_n);
  // fft_n elements is enough for both r2c output (needs N/2+1) and c2c output (needs N)
  chan->spectrum.fft_out = fftwf_alloc_complex(chan->spectrum.fft_n);
  // Wideband mode with real front end; use real->complex FFT
  if(chan->frontend->isreal)
    setup_real_fft(chan);
  else
    setup_complex_fft(chan);
}

static void setup_narrowband(struct channel *chan){
  // Set up narrow band (downconvert) mode
  double const blockrate = 1. / Blocktime; // Typically 50 Hz

  int const L = chan->frontend->L;
  int const M = chan->frontend->M;
  int const N = L + M - 1;

  double const margin = 400; // Allow 400 Hz for filter skirts at edge of I/Q receiver - calculate this
  unsigned long const samprate_base = lcm(lrint(blockrate),lrint(L*blockrate/N)); // Samprate must be allowed by receiver
  chan->spectrum.fft_n = lrint(chan->spectrum.bin_count + margin / chan->spectrum.rbw); // Minimum for search to avoid receiver filter skirt
  // This (int) cast should be cleaned up
  while(chan->spectrum.fft_n < 65536 && (!goodchoice(chan->spectrum.fft_n) || lrint(chan->spectrum.fft_n * chan->spectrum.rbw) % samprate_base != 0))
    chan->spectrum.fft_n++;
  chan->output.samprate = lrint(chan->spectrum.fft_n * chan->spectrum.rbw);
  chan->output.channels = 2; // IQ mode
  if(Verbose > 1)
    fprintf(stderr,"%s narrow spectrum: center %'.3lf Hz bin count %u, rbw %.1lf Hz, samprate %u Hz fft size %u\n",
	    chan->name,chan->tune.freq,chan->spectrum.bin_count,chan->spectrum.rbw,chan->output.samprate,chan->spectrum.fft_n);

  int blocklen = lrint(chan->output.samprate/blockrate);

  // Set up downconverter
  delete_filter_output(&chan->filter.out);
  int r = create_filter_output(&chan->filter.out,&chan->frontend->in,NULL,blocklen,COMPLEX);
  (void)r;
  assert(r == 0);

  chan->filter.max_IF = (double)(chan->output.samprate - margin)/2;
  chan->filter.min_IF = -chan->filter.max_IF;
  chan->filter2.blocking = 0; // Not used in this mode, make sure it's 0
  set_filter(&chan->filter.out,chan->filter.min_IF,chan->filter.max_IF,chan->filter.kaiser_beta);
  chan->filter.remainder = NAN; // Force init of downconverter
  chan->filter.bin_shift = 1010101010; // Unlikely - but a kludge, force init of phase rotator

  // Allocate persistent FFT I/O buffers sized for this fft_n; avoids per-call heap churn
  fftwf_free(chan->spectrum.fft_in_r);
  chan->spectrum.fft_in_r = NULL; // narrowband always uses complex FFT
  fftwf_free(chan->spectrum.fft_in_c);
  chan->spectrum.fft_in_c = fftwf_alloc_complex(chan->spectrum.fft_n);
  fftwf_free(chan->spectrum.fft_out);
  chan->spectrum.fft_out = fftwf_alloc_complex(chan->spectrum.fft_n);

  setup_complex_fft(chan);
}
static void setup_real_fft(struct channel *chan){
  if(chan->spectrum.plan != NULL)
    fftwf_destroy_plan(chan->spectrum.plan);
  float *in = fftwf_alloc_real(chan->spectrum.fft_n);
  assert(in != NULL);
  float complex *out = fftwf_alloc_complex(chan->spectrum.fft_n/2+1); // N/2 + 1 output points for real->complex
  assert(out != NULL);
  chan->spectrum.plan = plan_r2c(chan->spectrum.fft_n, in, out);
  fftwf_free(in);
  fftwf_free(out);
  assert(chan->spectrum.plan != NULL);
}
static void setup_complex_fft(struct channel *chan){
  // Wideband mode with complex front end, or narrowband mode with either front end
  if(chan->spectrum.plan != NULL)
    fftwf_destroy_plan(chan->spectrum.plan);

  float complex *in = fftwf_alloc_complex(chan->spectrum.fft_n);
  assert(in != NULL);
  float complex *out = fftwf_alloc_complex(chan->spectrum.fft_n);
  assert(out != NULL);
  chan->spectrum.plan = plan_complex(chan->spectrum.fft_n, in, out, FFTW_FORWARD);
  fftwf_free(in);
  fftwf_free(out);
  assert(chan->spectrum.plan != NULL);
}


#if RICE
// Experiment with rice coding of delta values in bin data
static void rice(struct channel *chan){
  assert(chan != NULL);

  // make vector of quantized measurements
  int const bin_count = chan->spectrum.bin_count;
  double scale = 1./chan->spectrum.step;
  int data[bin_count];
  int wbin = bin_count/2; // nyquist freq is most negative
  for(int i=0; i < bin_count; i++){
    double x = scale * (power2dB(chan->spectrum.bin_data[wbin++]) - chan->spectrum.base);
    if(x < 0)
      x = 0;
    data[i] = llrint(x);
    if(wbin == bin_count)
      wbin = 0;  // Continuing through dc and positive frequencies
  }

  int best_k = -1;
  int best_bits = 8 * bin_count;
  bool delta = false;

  // Test non-delta encoding
  for(int k=1;k < 6; k++){
    int bits = 0;
    int M = 1 << k;

    for(int i=0;i < bin_count; i++){
      int value = data[i]; // values are non-negative, no sign bit needed

      // Rice encode
      int q = value / M;
      int r = value % M;
      (void)r;

      // emit q 0-bits and a 1 bit
      bits += q + 1;
      // Emit r using k bits
      bits += k;
    }
    if(bits < best_bits){
      best_bits = bits;
      best_k = k;
    }
  }
  // Now test delta encoding
  for(int k=1;k < 6; k++){
    int bits = 0;
    int M = 1 << k;
    int prev = 0;

    for(int i=0;i < bin_count; i++){
      int delta = data[i] - prev;
      prev = data[i];
      int value = (abs(delta) << 1) | (delta < 0); // zig-zag encoding

      // Rice encode
      int q = value / M;
      int r = value % M;
      (void)r;

      // emit q 0-bits and a 1 bit
      bits += q + 1;
      // Emit r using k bits
      bits += k;
    }
    if(bits < best_bits){
      best_bits = bits;
      best_k = k;
      delta = true;
    }
  }
  fprintf(stderr,"rice encoding, delta %d,  k=%d: %d bits (%d bytes), %.1lf%%\n",delta,best_k,best_bits,best_bits/8,
	  100. * best_bits / (8 * bin_count));
}
#endif
