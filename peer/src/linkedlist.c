#include<stdlib.h>
#include<stdio.h>

#include "peer.h"

#define ITEM 7

int main() {
		peerlist_t *curr, *next, *head = NULL;
		int i;

		for(i=1;i<=10;i++) {
				curr = (peerlist_t *)malloc(sizeof(peerlist_t));
				curr->host.port = i;
				curr->next  = head;
				head = curr;
		}

		/* test the first item on the list */
		if (head && head->host.port==ITEM) {
				/* delete the first item in the list */
				curr = head->next;
				FREE(head);
				head = curr;
		}

		/* traverse the list */
		curr = head;
		while(curr) {
				/* test the next item on the list */
				next = curr->next;
				if (next && next->host.port==ITEM) {
						/* delete the next item in the list */
						curr->next = next->next;
						FREE(next);
						break;
				}
				curr = curr->next;
		}

		curr = head;
		while(curr) {
				printf("%d\n", curr->host.port);
				curr = curr->next ;
		}

}

/*
   struct list_el {
   int val;
   struct list_el * next;
   };

   typedef struct list_el item;

   void main() {
   item * curr, * head;
   int i;

   head = NULL;

   for(i=1;i<=10;i++) {
   curr = (item *)malloc(sizeof(item));
   curr->val = i;
   curr->next  = head;
   head = curr;
   }

   curr = head;

   while(curr) {
   printf("%d\n", curr->val);
   curr = curr->next ;
   }
   }
 */

