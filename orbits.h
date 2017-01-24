#define ERR_INVALID_INPUT 1
#define ERR_MALLOC_FAIL 2
#define ERR_FILE_OPEN 3
#define ERR_PGPLOT 4


// newtypes
typedef struct {
    float dt;
    float t;
    float G;
} environment;
typedef struct {
    float x;
    float y;
    float z;
} vector;
typedef struct {
    vector r;
    vector rold;
    vector v;
    vector a;
    float m;
} particle;

void orbitEulerStep(environment *env, particle *p1, particle *p2);
void orbitRKStep(environment *env, particle *p1, particle *p2);
void orbitLeapfrogStep(environment *env, particle *p1, particle *p2);

float dvxdtGrav(float t, float x, float y);
float dvydtGrav(float t, float x, float y);

void initParticleZeroes(particle *p);

// miscellaneous
void errorCase(const int errorCode);
