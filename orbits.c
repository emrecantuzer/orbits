#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <plplot/plplot.h>
#include "orbits.h"

int main(){
    printf("\n=============================================================\n");
    printf("This program is able to simulate a one-body problem orbiting\n");
    printf("the origin\n");
    printf("============================================================\n");

    //==========================================================================
    //-----------------------INITIALIZATIONS------------------------------------
    //==========================================================================


    //==========================================================================
    //-----------------------VARIABLE INITIALIZATIONS---------------------------
    //==========================================================================

    // generic counters
    int i, j;

    // physical system
    environment env;

    // particles
    particle pE, pRK, pLF, sun;

    //==========================================================================
    //-----------------------INITIALIZE VARIABLE VALUES-------------------------
    //==========================================================================
    env.t = 0.0;
    env.dt = 0.9;
    env.G = 1.0;

    initParticleZeroes(&pE);
    pE.r.x = 1.0;
                    pE.v.y = 1.0;
    initParticleZeroes(&pRK);
    pRK.r.x = 1.0;
                    pRK.v.y = 1.0;
    initParticleZeroes(&pLF);
    pLF.r.x = 1.0;
                    pLF.v.y = 1.0;
    initParticleZeroes(&sun);

    plsdev("xcairo");
    //plsetopt("drvopt","image_buffering=1");
    plinit();
    plenv(-5.0, 5.0, -5.0, 5.0, 0, 0);
    pllab("x", "y", "Title");

    //==========================================================================
    //-----------------------MAIN LOOP------------------------------------------
    //==========================================================================
    i=0;
    while(1) {
        orbitEulerStep(&env, &pE, &sun);
        orbitRKStep(&env, &pRK, &sun);
        orbitLeapfrogStep(&env, &pLF, &sun);

        //orbitRKStep(dt, &p1, &sun);

        env.t += env.dt;

        //printf("%f: %f, %f\n", env.t, pRK.r.x, pRK.r.y);

        //xs[i] = p1.r.x;
        //ys[i] = p1.r.y;

        plcol0(1);
        pljoin(pE.rold.x, pE.rold.y, pE.r.x, pE.r.y);
        plcol0(2);
        pljoin(pRK.rold.x, pRK.rold.y, pRK.r.x, pRK.r.y);
        plcol0(4);
        pljoin(pLF.rold.x, pLF.rold.y, pLF.r.x, pLF.r.y);
        plflush();

        env.t += env.dt;
        // slow down
        for (j=0; j<10000000; j++){}

        i++;
    }

    plend();
}


float dvxdtGrav(float t, float x, float y){
    float r = pow(x*x + y*y, 0.5);
    return (-1.0 * x) / (r*r*r);
}
float dvydtGrav(float t, float x, float y){
    float r = pow(x*x + y*y, 0.5);
    return (-1.0 * y) / (r*r*r);
}

void orbitLeapfrogStep(environment *env, particle *p1, particle *p2){
    float t = env->t; float dt = env->dt;

    // archive positions
    p1->rold.x = p1->r.x;
    p1->rold.y = p1->r.y;

    // vector distances
    float vecX = p1->r.x - p2->r.x;
    float vecY = p1->r.y - p2->r.y;
    // kick velocity forward half timestep
    p1->v.x += dvxdtGrav(t, vecX, vecY) * p2->m * env->G * dt/2;
    p1->v.y += dvydtGrav(t, vecX, vecY) * p2->m * env->G * dt/2;
    // drift position forward by full timestep
    p1->r.x += p1->v.x*dt;
    p1->r.y += p1->v.y*dt;
    // update vector distances
    vecX = p1->r.x - p2->r.x;
    vecY = p1->r.y - p2->r.y;
    // kick velocity forward by half timestep
    p1->v.x += dvxdtGrav(t, vecX, vecY) * p2->m * env->G * dt/2;
    p1->v.y += dvydtGrav(t, vecX, vecY) * p2->m * env->G * dt/2;
}

void orbitRKStep(environment *env, particle *p1, particle *p2){
    float t = env->t; float dt = env->dt;

    // init RK coefficients
    float kvx1; float kvy1;
    float kvx2; float kvy2;
    float kvx3; float kvy3;
    float kvx4; float kvy4;

    float krx1; float kry1;
    float krx2; float kry2;
    float krx3; float kry3;
    float krx4; float kry4;

    // vector distances
    float vecX = p1->r.x - p2->r.x;
    float vecY = p1->r.y - p2->r.y;

    // calculate forces (acceleration)
    krx1 = p1->v.x;
    kry1 = p1->v.y;
    kvx1 = dvxdtGrav(t, vecX, vecY) * p2->m * env->G;
    kvy1 = dvydtGrav(t, vecX, vecY) * p2->m * env->G;

    krx2 = p1->v.x + kvx1*dt/2.0;
    kry2 = p1->v.y + kvy1*dt/2.0;
    kvx2 = dvxdtGrav(t + dt/2.0, vecX + krx1*dt/2.0, vecY + kry1*dt/2.0) * p2->m * env->G;
    kvy2 = dvydtGrav(t + dt/2.0, vecX + krx1*dt/2.0, vecY + kry1*dt/2.0) * p2->m * env->G;

    krx3 = p1->v.x + kvx2*dt/2.0;
    kry3 = p1->v.y + kvy2*dt/2.0;
    kvx3 = dvxdtGrav(t + dt/2.0, vecX + krx2*dt/2.0, vecY + kry2*dt/2.0) * p2->m * env->G;
    kvy3 = dvydtGrav(t + dt/2.0, vecX + krx2*dt/2.0, vecY + kry2*dt/2.0) * p2->m * env->G;

    krx4 = p1->v.x + kvx3*dt;
    kry4 = p1->v.y + kvy3*dt;
    kvx4 = dvxdtGrav(t + dt, vecX + krx3*dt, vecY + kry3*dt) * p2->m * env->G;
    kvy4 = dvydtGrav(t + dt, vecX + krx3*dt, vecY + kry3*dt) * p2->m * env->G;

    //printf("k1: %f, %f\nk2: %f, %f\nk3: %f, %f\nk4: %f, %f\n", kvx1, kvy1, kvx2, kvy2, kvx3, kvy3, kvx4, kvy4);
    //printf("k1: %f, %f\nk2: %f, %f\nk3: %f, %f\nk4: %f, %f\n", krx1, kry1, krx2, kry2, krx3, kry3, krx4, kry4);

    // archive position
    p1->rold.x = p1->r.x;
    p1->rold.y = p1->r.y;

    // update velocity
    p1->v.x += ((kvx1 + 2.0*kvx2 + 2.0*kvx3 + kvx4) * dt/6.0);
    p1->v.y += ((kvy1 + 2.0*kvy2 + 2.0*kvy3 + kvy4) * dt/6.0);
    // update position
    p1->r.x += ((krx1 + 2.0*krx2 + 2.0*krx3 + krx4) * dt/6.0);
    p1->r.y += ((kry1 + 2.0*kry2 + 2.0*kry3 + kry4) * dt/6.0);
}


void orbitEulerStep(environment *env, particle *p1, particle *p2){
    float vecX = p1->r.x - p2->r.x;
    float vecY = p1->r.y - p2->r.y;
    // calculate forces (acceleration)
    p1->a.x = env->G * dvxdtGrav(env->t, vecX, vecY) * p2->m;
    p1->a.y = env->G * dvydtGrav(env->t, vecX, vecY) * p2->m;

    // archive position
    p1->rold.x = p1->r.x;
    p1->rold.y = p1->r.y;

    // update velocity
    p1->v.x += p1->a.x * env->dt;
    p1->v.y += p1->a.y * env->dt;
    // update position
    p1->r.x += p1->v.x * env->dt;
    p1->r.y += p1->v.y * env->dt;
}

void initParticleZeroes(particle *p){
    p->r.x = 0.0;    p->r.y = 0.0;    p->r.z = 0.0;      // position
    p->v.x = 0.0;    p->v.y = 0.0;    p->v.z = 0.0;      // velocity
    p->a.x = 0.0;    p->a.y = 0.0;    p->a.z = 0.0;      // acceleration

    p->rold.x = 0.0; p->rold.y = 0.0; p->rold.z = 0.0;
    p->m = 1.0;
}

void errorCase(const int errorCode){
    system("cat nagato");
    switch(errorCode){
        case ERR_INVALID_INPUT:
			printf("Error: invalid input\n");
			exit(-1);
		case ERR_MALLOC_FAIL:
			printf("Error: out of memory\n");
			exit(-1);
		case ERR_FILE_OPEN:
			printf("Error: file cannot be opened\n");
			exit(-1);
		case ERR_PGPLOT:
			printf("Error: cannot open pgplot window\n");
			exit(-1);
    }
}
