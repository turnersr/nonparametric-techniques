/* foo.h */

//#include <cstddef>

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  void RFFT(double* x, int array_size, double* return_y);

  void RIFFT(double* x, int array_size, double* return_y);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */
