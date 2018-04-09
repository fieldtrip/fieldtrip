#include <limits.h>

#define ACQ_MSGQ_SIZE      600
#define ACQ_MSGQ_SHMKEY    0x39457f73
#define ACQ_MSGQ_SHMPROJID 12345
#define ACQ_MSGQ_SHMPATH   "/opt/ctf/bin/Acq"
#define ACQ_BUFFER_SIZE    40000

/*
 * The ACQ_BUFFER_SIZE should be set to 28160 for older acquisition software and to 40000
 * for newer (beta) versions of the acquisition software from approximately 2016 onwards.
 * See also http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3185
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

