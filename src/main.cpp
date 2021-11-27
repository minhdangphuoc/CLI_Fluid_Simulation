#include <iostream>
#include <cstdlib>
#include <math.h>

int N = 16;
int iter = 4;

#define IX(x, y) ((x) + (y) * N)

struct FluidCube {
    int size;
    float iter;
    float diff;
    float visc;
    
    float *s;
    float *density;
    
    float *Vx;
    float *Vy;

    float *Vx0;
    float *Vy0;
};


typedef struct FluidCube FluidCube;

FluidCube *FluidCubeCreate(int size, int diffusion, int viscosity, float iter)
{
    FluidCube *cube = new FluidCube;
    int N = size;
    
    cube->size = size;
    cube->iter = iter;
    cube->diff = diffusion;
    cube->visc = viscosity;
    
    cube->s = (float*)calloc(N * N * N, sizeof(float));
    cube->density = (float*)calloc(N * N * N, sizeof(float));
    
    cube->Vx = (float*)calloc(N * N * N, sizeof(float));
    cube->Vy = (float*)calloc(N * N * N, sizeof(float));
    
    cube->Vx0 = (float*)calloc(N * N * N, sizeof(float));
    cube->Vy0 = (float*)calloc(N * N * N, sizeof(float));
    
    return cube;
}

void FluidCubeFree(FluidCube *cube)
{
    free(cube->s);
    free(cube->density);
    
    free(cube->Vx);
    free(cube->Vy);

    free(cube->Vx0);
    free(cube->Vy0);
    
    free(cube);
}


static void set_bnd(int b, float *x)
{
    for (int i = 1; i < N - 1; i++) {
        x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N - 1)] = b == 2 ? -x[IX(i, N - 2)] : x[IX(i, N - 2)];
    }
    for (int j = 1; j < N - 1; j++) {
        x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
        x[IX(N - 1, j)] = b == 1 ? -x[IX(N - 2, j)] : x[IX(N - 2, j)];
    }
    
    x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N - 1)] = 0.5 * (x[IX(1, N - 1)] + x[IX(0, N - 2)]);
    x[IX(N - 1, 0)] = 0.5 * (x[IX(N - 2, 0)] + x[IX(N - 1, 1)]);
    x[IX(N - 1, N - 1)] = 0.5 * (x[IX(N - 2, N - 1)] + x[IX(N - 1, N - 2)]);
}

static void lin_solve(int b, float *x, float *x0, float a, float c)
{
    float cRecip = 1.0 / c;
    for (int k = 0; k < iter; k++) {
            for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                    x[IX(i, j)] =
                    (x0[IX(i, j)] +
                        a *
                        (x[IX(i + 1, j)] +
                            x[IX(i - 1, j)] +
                            x[IX(i, j + 1)] +
                            x[IX(i, j - 1)])) * cRecip;
                }
            }
        }
        set_bnd(b, x);
}

static void diffuse (int b, float *x, float *x0, float diff)
{
    float a = iter * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a);
}

static void advect(int b, float *d, float *d0,  float *velocX, float *velocY)
{
    float i0, i1, j0, j1;
    
    float dtx = iter * (N - 2);
    float dty = iter * (N - 2);
    
    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;
    
    float Nfloat = N - 2;
    float ifloat, jfloat;
    int i, j;
    
        for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
            for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                tmp1 = dtx * velocX[IX(i, j)];
                tmp2 = dty * velocY[IX(i, j)];
                x    = ifloat - tmp1; 
                y    = jfloat - tmp2;
                
                if(x < 0.5f) x = 0.5f; 
                if(x > Nfloat + 0.5f) x = Nfloat + 0.5f; 
                i0 = floorf(x); 
                i1 = i0 + 1.0f;
                if(y < 0.5f) y = 0.5f; 
                if(y > Nfloat + 0.5f) y = Nfloat + 0.5f; 
                j0 = floorf(y);
                j1 = j0 + 1.0f; 
                
                s1 = x - i0; 
                s0 = 1.0f - s1; 
                t1 = y - j0; 
                t0 = 1.0f - t1;
                
                int i0i = i0;
                int i1i = i1;
                int j0i = j0;
                int j1i = j1;

                
                d[IX(i, j)] =
                    s0 * (t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)]) +
                    s1 * (t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)]);
            }
        }
    
    set_bnd(b, d);
}

static void project(float *velocX, float *velocY, float *p, float *div)
{
    for (int k = 1; k < N - 1; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                div[IX(i, j)] = -0.5f*(
                         velocX[IX(i+1, j)]
                        -velocX[IX(i-1, j)]
                        +velocY[IX(i  , j+1)]
                        -velocY[IX(i  , j-1)]
                    )/N;
                p[IX(i, j)] = 0;
            }
        }
    }
    set_bnd(0, div); 
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 6);

        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                velocX[IX(i, j)] -= 0.5f * (  p[IX(i+1, j)]
                                                -p[IX(i-1, j)]) * N;
                velocY[IX(i, j)] -= 0.5f * (  p[IX(i, j+1)]
                                                -p[IX(i, j-1)]) * N;
            }
        }
    set_bnd(1, velocX);
    set_bnd(2, velocY);
}

void FluidCubeStep(FluidCube *cube) // in Main
{
    int N          = cube->size;
    float visc     = cube->visc;
    float diff     = cube->diff;
    float *Vx      = cube->Vx;
    float *Vy      = cube->Vy;
    float *Vx0     = cube->Vx0;
    float *Vy0     = cube->Vy0;
    float *s       = cube->s;
    float *density = cube->density;
    
    diffuse(1, Vx0, Vx, visc);
    diffuse(2, Vy0, Vy, visc);

    project(Vx0, Vy0, Vx, Vy);
    
    advect(1, Vx, Vx0, Vx0, Vy0);
    advect(2, Vy, Vy0, Vx0, Vy0);
    
    project(Vx, Vy, Vx0, Vy0);
    
    diffuse(0, s, density, diff);
    advect(0, density, s, Vx, Vy);
}

void FluidCubeAddDensity(FluidCube *cube, int x, int y, float amount)
{
    int N = cube->size;
    cube->density[IX(x, y)] += amount;
}

void FluidCubeAddVelocity(FluidCube *cube, int x, int y, float amountX, float amountY)
{
    int N = cube->size;
    int index = IX(x, y);
    
    cube->Vx[index] += amountX;
    cube->Vy[index] += amountY;
}

void printMap(FluidCube *cube)
{
    float d;
    char map[N][N];
    for (int i = 0; i < N+2; i++){
        std :: cout << '/';
    }
    std :: cout << std :: endl;
    for (int i = 0; i < 16; i++)
    {
        std :: cout << '/';
        for (int j = 0; j < 16; j++)
        {
            d = cube->density[IX(i, j)];
            if (d > 7) std :: cout << "▓";
            else if (d > 3) std :: cout << "▒";
            else if (d > 0) std :: cout << "░";
            else std :: cout << ' ';
        }
        std :: cout << '/';
        std :: cout << std :: endl;
    }
    for (int i = 0; i < N+2; i++){
        std :: cout << '/';
    }
}

void UIPrint(int size, FluidCube *cube)
{
    std::cout << "Refresh" << std::endl;
    float x = 7;
    float y = 7;
    FluidCubeAddDensity(cube, x, y, 100);
    FluidCubeAddVelocity(cube, x, y, x - rand() % 16, y - rand() % 16);
    printMap(cube);
    std::cin.get();
    system("clear");
}

int main() 
{
    FluidCube *cube = FluidCubeCreate(16, 0.2, 0, 0.0000001);
    while (true)
    {
        FluidCubeStep(cube);
        UIPrint(16, cube);
    }
        
}