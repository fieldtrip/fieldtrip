/* binsgeq
 *
 * binary search of a sorted (in ascending order) double array <*vec>
 * of length <len> for a first element which is greater than or equal
 * to a key <item>. If no greater element is found, returns <len>.
 */


/* Copyright (c) 1998-2004 Aki Vehtari  
 *
 *This software is distributed under the GNU General Public 
 *License (version 3 or later); please refer to the file 
 *License.txt, included with the software, for details.
 *
 */

double binsgeq(double *vec, int len, double item)
{
  int low=0, high=len-1, mid;
  if (item<=vec[0]) {
    return 1.0; 
  }
  while ((high-low) > 1) {
    mid=(high+low)/2;
    if (item>vec[mid])
      low = mid;
    else
      high = mid;
  }
  return high+1.0; 
}
