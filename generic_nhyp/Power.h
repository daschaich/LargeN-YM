#define P4(x)    (x)*(x)*(x)*(x)
#define Power(x,n) (n==1 ? (x) :                \
                   (n==2 ? (x)*(x) :            \
                   (n==3 ? (x)*(x)*(x) :        \
                   (n==4 ? P4(x) :              \
                   (n==5 ? P4(x)*(x) :          \
                   (n==6 ? P4(x)*(x)*(x) :      \
                   (n==7 ? P4(x)*(x)*(x)*(x) :  \
                   (n==8 ? P4(x)*P4(x) :        \
                   (n==9 ? P4(x)*P4(x)*(x) : -1 )))))))))
