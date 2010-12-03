#ifndef _dacq_message_h
#define _dacq_message_h

typedef struct {
  int     kind;			/* What is this data? */
  int     type;			/* What is its type */
  int     size;			/* Size of item */
  int     loc;			/* Position in file */
  int     shmem_buf;		/* Shared mem block */
  int     shmem_loc;		/* Not used, set to -1 */
} dacqDataMessageRec,*dacqDataMessage;

#define DATA_MESS_SIZE sizeof(dacqDataMessageRec)

#endif
