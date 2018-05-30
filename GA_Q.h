struct particle
{
    double x1,x2,x3;
    double fitness;
    double frecuencia;
    double amplitud;
    double v1,v2,v3;
    double best_x1,best_x2,best_x3;   
    double best_p_f;
};

struct Best_Particle
{
    double x1,x2,x3;
    double fitness;
    double frecuencia;
    double amplitud;
    int generation;
};

struct Population
{   
    double x1,x2,x3;
    double fitness;
    double frecuencia;
    double amplitud;
};

///////////////////////////////* Functions prototypes *///////////////////////////
void Allocate_Memory(int popsize);
void nomemory(char *string);
void Free_All();
void Get_Previous_Individuals();
void Initialize_Population(int popsize);
void Evaluate_Individual(struct particle *ind);
void randomize(float seed);
void warmup_random(float random_seed);
void advance_random();
float rndreal(float lo ,float hi);
float randomperc();
void add_individual(struct particle * ind);
int buscar (struct particle * ind);
void show_pop(int popsize);
void Initial_Report(int popsize,int gmax,float personal,float global,float i_max,float i_min,float seed);
void Update_Bp_Gg(int popsize,int gen);
void Update_Position(int popsize,float p_pond,float g_pond,float inertia);
double round_number(double number);
void statistics(struct particle *pop,int gen,int popsize,int variables);
void report(float min,int peor,float max,int mejor,float mean_gen,int gen,float ecm);
void Track_Particles_Function(struct particle *pop,int popsize);
double ECM(struct particle *pop,int popsize,int variables);

//////////////////////////////////////////////////////////////////////////////////
struct particle *initial_pop;
struct particle *oldpop;                    /* last generation of individuals */
struct Population *globalpop;                    
struct Best_Particle bestp;
struct Best_Particle bestp_neighborhood;
double sumfitness;                      /* summed fitness for entire population */
double oldrand[55];                               /* Array of 55 random numbers */
int jrand;                                             /* current random number */
double rndx2;                                /* used with random normal deviate */
int rndcalcflag;                             /* used with random normal deviate */
