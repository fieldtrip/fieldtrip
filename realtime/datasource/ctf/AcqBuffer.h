#include <limits.h>

#define ACQ_MSGQ_SIZE      600
#define ACQ_MSGQ_SHMKEY    0x39457f73
#define ACQ_MSGQ_SHMPROJID 12345
#define ACQ_MSGQ_SHMPATH   "/opt/ctf/bin/Acq"

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
  int data[28160];
} ACQ_MessagePacketType;

