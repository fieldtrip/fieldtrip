#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "message.h"

int main(void) {
    printf("sizeof (CHAR        ) = %d, expected %d\n", WORDSIZE_CHAR        ,1);
    printf("sizeof (UINT8       ) = %d, expected %d\n", WORDSIZE_UINT8       ,1);
    printf("sizeof (UINT16      ) = %d, expected %d\n", WORDSIZE_UINT16      ,2);
    printf("sizeof (UINT32      ) = %d, expected %d\n", WORDSIZE_UINT32      ,4);
    printf("sizeof (UINT64      ) = %d, expected %d\n", WORDSIZE_UINT64      ,8);
    printf("sizeof (INT8        ) = %d, expected %d\n", WORDSIZE_INT8        ,1);
    printf("sizeof (INT16       ) = %d, expected %d\n", WORDSIZE_INT16       ,2);
    printf("sizeof (INT32       ) = %d, expected %d\n", WORDSIZE_INT32       ,4);
    printf("sizeof (INT64       ) = %d, expected %d\n", WORDSIZE_INT64       ,8);
    printf("sizeof (FLOAT32     ) = %d, expected %d\n", WORDSIZE_FLOAT32     ,4);
    printf("sizeof (FLOAT64     ) = %d, expected %d\n", WORDSIZE_FLOAT64     ,8);
    printf("\n");
    printf("sizeof (message_t   ) = %d, expected %d\n", sizeof(message_t)    ,8);
    printf("sizeof (header_t    ) = %d, expected %d\n", sizeof(header_t)     ,8);
    printf("sizeof (data_t      ) = %d, expected %d\n", sizeof(data_t)       ,8);
    printf("sizeof (event_t     ) = %d, expected %d\n", sizeof(event_t)      ,8);
    printf("sizeof (datasel_t   ) = %d, expected %d\n", sizeof(datasel_t)    ,8);
    printf("sizeof (eventsel_t  ) = %d, expected %d\n", sizeof(eventsel_t)   ,8);
    printf("\n");
    printf("sizeof (messagedef_t) = %d, expected %d\n", sizeof(messagedef_t) ,8);
    printf("sizeof (headerdef_t ) = %d, expected %d\n", sizeof(headerdef_t)  ,24);
    printf("sizeof (datadef_t   ) = %d, expected %d\n", sizeof(datadef_t)    ,16);
    printf("sizeof (eventdef_t  ) = %d, expected %d\n", sizeof(eventdef_t)   ,32);
}

