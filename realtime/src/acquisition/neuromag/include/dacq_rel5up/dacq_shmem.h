#ifndef _dacq_shmem_h
#define _dacq_shmem_h

#define SHM_FILE       "/neuro/dacq/shmem/data_server"
#define SHM_FAIL_FILE  "/neuro/dacq/raw/data_server_shmem"
#define SHM_MAX_CLIENT 10

#define SHM_MAX_DATA   500*1500*4
#define SHM_NUM_BLOCKS 100
#define SHM_NO_BUF     -1

typedef struct {
  int client_id;
  int done;
} *dacqShmClient,dacqShmClientRec;

typedef struct {
  dacqShmClientRec clients[SHM_MAX_CLIENT];
  unsigned char data[SHM_MAX_DATA];
} *dacqShmBlock,dacqShmBlockRec;

#define SHM_SIZE SHM_NUM_BLOCKS*sizeof(dacqShmBlockRec)

#endif

