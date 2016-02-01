#include <siemensap.h>
#include <stdio.h>

char buffer[1024*1024];

int main(int argc, char *argv[])
{
   FILE *fp;
   size_t len;
   sap_item_t *S;
   sap_essentials_t SE;
   
   if (argc < 2) 
   {
      printf("Usage: tsap <mrprot.txt>\n");
      return 1;
   }
   
   fp = fopen(argv[1], "rb");
   if (fp == NULL)
   {
      printf("Files '%s' cannot be opened.\n", argv[1]);
      return 1;
   }
   
   len = fread(buffer, 1, sizeof(buffer), fp);
   
   printf("File contains %i bytes\n", len);
   fclose(fp);
   
   S = sap_parse(buffer, len);
   
   printf("\n\nReversing\n\n");
   
   sap_reverse_in_place(&S);
   
   printf("\n\n\n\n");
   sap_print(S);
   
   sap_get_essentials(S, &SE);
   
   printf("TR = %i\n", SE.TR);
   printf("Res: %i x %i x %i\n", SE.readoutPixels, SE.phasePixels, SE.numberOfSlices);
   printf("Contrasts: %i\n", SE.numberOfContrasts);
   printf("FOV: %f x %f\n", SE.readoutFOV, SE.phaseFOV);
   printf("Thickn: %f\n", SE.sliceThickness);
	
   sap_destroy(S);
}