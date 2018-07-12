#include <limits.h>

#define ACQ_MSGQ_SIZE      600
#define ACQ_MSGQ_SHMKEY    0x39457f73
#define ACQ_MSGQ_SHMPROJID 12345
#define ACQ_MSGQ_SHMPATH   "/opt/ctf/bin/Acq"

#define ACQ_BUFFER_SIZE    28160 // this is needed for all other software versions
// #define ACQ_BUFFER_SIZE 40000 // this is needed for software version 6.1.5-el6_7.x86_64-20160720-3344 

/*
 * See also 
 *   http://www.fieldtriptoolbox.org/development/realtime/ctf#different_software_versions
 *   http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3185
 *   https://github.com/fieldtrip/fieldtrip/issues/699
 *   https://github.com/fieldtrip/fieldtrip/issues/724
 */

typedef enum {
  ACQ_MSGQ_SETUP_COLLECTION,
  ACQ_MSGQ_DATA,
  ACQ_MSGQ_CLOSE_CONNECTION,
  ACQ_MSGQ_INVALID = INT_MAX
} ACQ_MessageType;

typedef struct {
  ACQ_MessageType message_type;
  int messageId;
  int sampleNumber;
  int numSamples;
  int numChannels;
  int data[ACQ_BUFFER_SIZE];
} ACQ_MessagePacketType;

