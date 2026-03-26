// read FFT bin energies from spectrum pseudo-demod and format similar to rtl_power - out of date
// Copyright 2023 Phil Karn, KA9Q

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#if defined(linux)
#include <bsd/string.h>
#endif
#include <assert.h>
#include <getopt.h>
#include <sysexits.h>
#include <fcntl.h>

#include "misc.h"
#include "status.h"
#include "multicast.h"
#include "radio.h"

struct sockaddr Metadata_dest_socket;      // Dest of metadata (typically multicast)
struct sockaddr Metadata_source_socket;      // Source of metadata
int IP_tos;
int Mcast_ttl = 1;
const char *App_path;
const char *Target;
int Verbose;
uint32_t Ssrc;
char Iface[1024]; // Multicast interface to talk to front end
int Status_fd = -1;
int64_t Timeout = BILLION; // Retransmission timeout
bool details;   // Output bin, frequency, power, newline
char const *Source;
struct sockaddr_storage *Source_socket;

static char const Optstring[] = "b:c:C:df:hi:o:s:t:T:vw:V";
static struct  option Options[] = {
  {"bins", required_argument, NULL, 'b'},
  {"count", required_argument, NULL, 'c'},
  {"details", no_argument, NULL, 'd'},
  {"frequency", required_argument, NULL, 'f'},
  {"help", no_argument, NULL, 'h'},
  {"interval", required_argument, NULL, 'i'},
  {"ssrc", required_argument, NULL, 's'},
  {"timeout", required_argument, NULL, 'T'},
  {"verbose", no_argument, NULL, 'v'},
  {"version", no_argument, NULL, 'V'},
  {"bin-width", required_argument, NULL, 'w'},
  {"crossover", required_argument, NULL, 'C'},
  {"source", required_argument, NULL, 'o'},
  {NULL, 0, NULL, 0},
};


int extract_powers(float *power,int npower,uint64_t *time,double *freq,double *rbw,int32_t const ssrc,uint8_t const * const buffer,size_t length);

void help(){
  fprintf(stderr,"Usage: %s [-v|--verbose] [-V|--version] [-f|--frequency freq] [-w|--bin-width rbw] [-b|--bins bins] [-c|--count count] [-i|--interval interval] [-T|--timeout timeout] [-d|--details] -s|--ssrc ssrc mcast_addr [-o|--source <source name-or-address>\n",App_path);
  exit(1);
}

int main(int argc,char *argv[]){
  App_path = argv[0];
  int count = 1;     // Number of updates. -1 means infinite
  double interval = 5; // Period between updates, sec
  double frequency = -1;
  int bins = 0;
  double rbw = 0;
  double crossover = -1;
  {
    int c;
    while((c = getopt_long(argc,argv,Optstring,Options,NULL)) != -1){
      switch(c){
      case 'b':
	bins = atoi(optarg);
	break;
      case 'c':
	count = atoi(optarg);
	break;
      case 'C':
	crossover = strtod(optarg,NULL);
	break;
      case 'd':
	details = true;
	break;
      case 'f':
	frequency = parse_frequency(optarg,true);
	break;
      case 'h':
	help();
	break;
      case 'i':
	interval = strtod(optarg,NULL);
	break;
      case 's':
	Ssrc = atoi(optarg); // Send to specific SSRC
	break;
      case 'T':
	Timeout = (int64_t)(BILLION * strtod(optarg,NULL)); // Retransmission timeout
	break;
      case 'v':
	Verbose++;
	break;
      case 'w':
	rbw = strtod(optarg,NULL);
	break;
      case 'V':
	VERSION();
	exit(EX_OK);
      case 'o':
	Source = optarg;
	break;
      default:
	fprintf(stdout,"Unknown option %c\n",c);
	help();
	break;
      }
    }
  }
  if(argc <= optind)
    help();

  Target = argv[optind];
  resolve_mcast(Target,&Metadata_dest_socket,DEFAULT_STAT_PORT,Iface,sizeof(Iface),0);
  if(Verbose)
    fprintf(stderr,"Resolved %s -> %s\n",Target,formatsock(&Metadata_dest_socket,false));

  if(Source != NULL){
    Source_socket = calloc(1,sizeof(struct sockaddr_storage));
    if(Verbose)
      fprintf(stdout,"Resolving source %s\n",Source);
    resolve_mcast(Source,Source_socket,0,NULL,0,0);
  }

  Status_fd = listen_mcast(Source_socket,&Metadata_dest_socket,Iface);
  if(Status_fd == -1){
    fprintf(stderr,"Can't listen to mcast status %s\n",Target);
    exit(1);
  }
  int Ctl_fd = output_mcast(&Metadata_dest_socket,Iface,Mcast_ttl,IP_tos);
  if(Ctl_fd == -1){
    fprintf(stderr,"connect to mcast control failed: %s\n",strerror(errno));
    exit(1);
  }

  // Send command to set up the channel?? Or do in a separate command? We'd like to reuse the same demod & ssrc,
  // which is hard to do in one command, as we'd have to stash the ssrc somewhere.
  while(true){
    uint8_t buffer[PKTSIZE];

    // Accumulate power over the interval by polling repeatedly
    double accum[PKTSIZE / sizeof(float)]; // double for accumulation precision
    memset(accum, 0, sizeof(accum));
    int polls = 0;
    size_t npower = 0;
    uint64_t time = 0;
    double r_freq = 0;
    double r_rbw = 0;
    int64_t const interval_end = gps_time_ns() + (int64_t)(interval * BILLION);

    do {
      // Build command with fresh tag each poll
      uint8_t *bp = buffer;
      *bp++ = 1; // Command
      encode_int(&bp,OUTPUT_SSRC,Ssrc);
      uint32_t tag = (uint32_t)random();
      encode_int(&bp,COMMAND_TAG,tag);
      encode_int(&bp,DEMOD_TYPE,SPECT_DEMOD);
      if(frequency >= 0)
	encode_double(&bp,RADIO_FREQUENCY,frequency);
      if(bins > 0)
	encode_int(&bp,BIN_COUNT,bins);
      if(rbw > 0)
	encode_float(&bp,RESOLUTION_BW,rbw);
      if(crossover >= 0)
	encode_float(&bp,CROSSOVER,crossover);
      // Set server-side integration count
      // With 50% overlap (default), each FFT step advances fft_n/2 samples = 1/(2*rbw) sec
      // So frames needed = interval / (1/(2*rbw)) = interval * rbw * 2
      if(rbw > 0 && interval > 0){
	int avg = (int)(interval * rbw * 2);
	if(avg < 1)
	  avg = 1;
	encode_int(&bp,SPECTRUM_AVG,avg);
      }
      encode_eol(&bp);
      ssize_t const cmd_len = bp - buffer;

      if(Verbose > 1){
	fprintf(stderr,"Sent:");
	dump_metadata(stderr,buffer+1,cmd_len-1,details ? true : false);
      }
      if(sendto(Ctl_fd, buffer, cmd_len, 0, &Metadata_dest_socket, sizeof Metadata_dest_socket) != cmd_len){
	perror("command send");
	usleep(10000);
	continue;
      }
      // Deadline must cover server-side integration time
      int64_t poll_timeout = Timeout;
      if(rbw > 0 && interval > 0){
	int64_t integration_ns = (int64_t)(interval * BILLION) + BILLION;
	if(integration_ns > poll_timeout)
	  poll_timeout = integration_ns;
      }
      int64_t deadline = gps_time_ns() + poll_timeout;
      ssize_t length = 0;
      do {
	fd_set fdset;
	FD_ZERO(&fdset);
	FD_SET(Status_fd,&fdset);
	int n = Status_fd + 1;
	int64_t timeout = deadline - gps_time_ns();
	if(timeout < 0)
	  timeout = 0;
	struct timespec ts;
	ns2ts(&ts,timeout);
	n = pselect(n,&fdset,NULL,NULL,&ts,NULL);
	if(n <= 0 && timeout == 0)
	  break; // timed out waiting for this poll
	if(!FD_ISSET(Status_fd,&fdset))
	  continue;
	socklen_t ssize = sizeof(Metadata_source_socket);
	length = recvfrom(Status_fd,buffer,sizeof(buffer),0,(struct sockaddr *)&Metadata_source_socket,&ssize);
      } while(length < 2 || (enum pkt_type)buffer[0] != STATUS || Ssrc != get_ssrc(buffer+1,length-1) || tag != get_tag(buffer+1,length-1));

      if(length < 2)
	continue; // timed out, try again

      if(Verbose > 1){
	fprintf(stderr,"Received:");
	dump_metadata(stderr,buffer+1,length-1,details ? true : false);
      }
      float poll_powers[PKTSIZE / sizeof(float)];
      size_t np = extract_powers(poll_powers,sizeof(poll_powers) / sizeof(poll_powers[0]),
				 &time,&r_freq,&r_rbw,Ssrc,buffer+1,length-1);
      if(np <= 0){
	if(Verbose)
	  fprintf(stderr,"Invalid response, length %lu\n",np);
	continue;
      }
      if(npower == 0)
	npower = np; // first valid response sets the bin count
      else if(np != npower)
	continue; // bin count changed, skip

      for(size_t i = 0; i < npower; i++)
	accum[i] += poll_powers[i];
      polls++;
    } while(gps_time_ns() < interval_end);

    if(polls == 0 || npower == 0){
      fprintf(stderr,"No valid responses during interval\n");
      goto again;
    }
    // Average the accumulated powers
    for(size_t i = 0; i < npower; i++)
      accum[i] /= polls;

    if(Verbose)
      fprintf(stderr,"Integrated %d polls over %.1f sec\n",polls,interval);

    // Note from VK5QI:
    // the output format from that utility matches that produced by rtl_power, which is:
    //2022-04-02, 16:24:55, 400050181, 401524819, 450.13, 296, -52.95, -53.27, -53.26, -53.24, -53.40, <many more points here>
    // date, time, start_frequency, stop_frequency, bin_size_hz, number_bins, data0, data1, data2

    char gps[1024];
    printf("%s,",format_gpstime_iso8601(gps,sizeof(gps),time));

    size_t const first_neg_bin = (npower + 1)/2;
    double base = r_freq - r_rbw * (npower/2);
    printf(" %.0lf, %.0lf, %.0lf, %lu",
	   base, base + r_rbw * (npower-1), r_rbw, npower);

    // Find lowest non-zero entry, use the same for zero power to avoid -infinity dB
    double lowest = INFINITY;
    for(size_t i=0; i < npower; i++){
      if(accum[i] < 0){
	fprintf(stderr,"Invalid power %g in response\n",accum[i]);
	goto again;
      }
      if(accum[i] > 0 && accum[i] < lowest)
	lowest = accum[i];
    }
    double const min_db = lowest != INFINITY ? power2dB(lowest) : 0;

    if (details){
      printf("\n");
      for(size_t i=first_neg_bin ; i < npower; i++){
        printf("%lu %lf %.2lf\n",i,base,(accum[i] == 0) ? min_db : power2dB(accum[i]));
        base += r_rbw;
      }
      for(size_t i=0; i < first_neg_bin; i++){
        printf("%lu %lf %.2lf\n",i,base,(accum[i] == 0) ? min_db : power2dB(accum[i]));
        base += r_rbw;
      }
    } else {
      for(size_t i= first_neg_bin; i < npower; i++)
        printf(", %.2lf",(accum[i] == 0) ? min_db : power2dB(accum[i]));
      for(size_t i=0; i < first_neg_bin; i++)
        printf(", %.2lf",(accum[i] == 0) ? min_db : power2dB(accum[i]));
    }
    printf("\n");
    if(--count == 0)
      break;
  again:;
  }
  exit(0);
}

// Decode only those status fields relevant to spectrum measurement
// Return number of bins
int extract_powers(float *power,int npower,uint64_t *time,double *freq,double *rbw,int32_t const ssrc,uint8_t const * const buffer,size_t length){
#if 0  // use later
  double l_lo1 = 0,l_lo2 = 0;
#endif
  int l_ccount = 0;
  uint8_t const *cp = buffer;
  int l_count = 0;

  while(cp < &buffer[length]){
    enum status_type const type = *cp++; // increment cp to length field

    if(type == EOL)
      break; // End of list

    unsigned int optlen = *cp++;
    if(optlen & 0x80){
      // length is >= 128 bytes; fetch actual length from next N bytes, where N is low 7 bits of optlen
      int length_of_length = optlen & 0x7f;
      optlen = 0;
      while(length_of_length > 0){
	optlen <<= 8;
	optlen |= *cp++;
	length_of_length--;
      }
    }
    if(cp + optlen >= buffer + length)
      break; // Invalid length
    switch(type){
    case EOL: // Shouldn't get here
      goto done;
    case GPS_TIME:
      *time = decode_int64(cp,optlen);
      break;
    case OUTPUT_SSRC: // Don't really need this, it's already been checked
      if((int32_t)decode_int32(cp,optlen) != ssrc)
	return -1; // Not what we want
      break;
    case DEMOD_TYPE:
      {
	const int i = decode_int(cp,optlen);
	if(i != SPECT_DEMOD)
	  return -1; // Not what we want
      }
      break;
    case RADIO_FREQUENCY:
      *freq = decode_double(cp,optlen);
      break;
#if 0  // Use this to fine-tweak freq later
    case FIRST_LO_FREQUENCY:
      l_lo1 = decode_double(cp,optlen);
      break;
    case SECOND_LO_FREQUENCY: // ditto
      l_lo2 = decode_double(cp,optlen);
      break;
#endif
    case BIN_DATA:
      l_count = optlen/sizeof(float);
      if(l_count > npower)
	return -2; // Not enough room in caller's array
      // Note these are still in FFT order
      for(int i=0; i < l_count; i++){
	power[i] = decode_float(cp + i * sizeof(float),sizeof(float));
      }
      break;
    case RESOLUTION_BW:
      *rbw = decode_float(cp,optlen);
      break;
    case BIN_COUNT: // Do we check that this equals the length of the BIN_DATA tlv?
      l_ccount = decode_int(cp,optlen);
      break;
    default:
      break;
    }
    cp += optlen;
  }
 done:
  ;
  if(l_ccount == 0 || l_count != l_ccount)
    return 0;
  return l_count;
}
