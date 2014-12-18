//
//**************************** General Utility Library Definitions ***********************
//

class Pair<F,S>
{
  F first;
  S second;
  
  Pair(F first, S second) {
    this.first = first;
    this.second = second;
  }
}
  
<T> void swap (Pair<T,T> pair) {
  T tmp = pair.second;
  pair.second = pair.first;
  pair.first = tmp;
}

<T extends Comparable<T>> T clamp (T value, T _min, T _max) {
  if (value.compareTo(_min) < 0) {
    return _min;
  } 
  if (value.compareTo(_max) > 0) {
    return _max;
  }  
  return value;
}

int pmod(int value, int mod) {
  int r = value % mod;   
  if (r < 0) r += mod; // Correct positive modulus
  return r;
}
