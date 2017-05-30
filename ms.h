struct devent {
	double time ;
	int popi ; // popID 1
	int popj ; // popID 2
	double paramv ; // parametet for event
	double **mat ;
	char detype ; // type of event, n, g, s, j or whatever
	struct devent *nextde ;
	} ;
	
struct c_params {
	int npop ; // # of subpopulation, initial: 1
	int nsam ; // # of samples
	int *config; // # of samples in each pop, initial: nsam
	double **mig_mat; // migration matrix, initial: 1x1 = 0
	double r ; // 4Nr, initial: 0
	int nsites ; // length , initial: 2
	double f;    // g/r, initial: 0
	double track_len ; // 1/q, initial: 0
	double *size; // population size, initial: 1
	double *alphag;
	struct devent *deventlist ; // demographic events, -e
    
    // N E W ========
    double f2 ; //g', initial: 0
    double track_len2 ; // length on gene convergion between species, initial: 0
    double divergence ; //divergence rate that can recombine with host genome, initial: 0, must be < 1.0
    // N E W =====end
	
	int n_two_dim ; //st no. of elements in two_dim
    int n_tractl ; //st no. of elements in tractl (tract length for cleft2)
    
    // N E W 2 ====== Size of 2-dim matrix (dist.txt)
    int size_tract ;
    int size_div ;
    // N E W 2 ===end
    
} ;
	
struct m_params {
	double theta ;  // 4Nu, initial: 0
	int segsitesin; // Number of seg sites, initial: 0
	int treeflag;   // -T tree output, initial: 0
	int timeflag;   // -L length of tree, initial: 0
	int mfreq ; // minor allele? initial: 1
} ;
	 
struct params { 
	struct c_params cp ;
	struct m_params mp ;
	int commandlineseedflag ;
} ;
	
void ranvec(int, double[]);// modified 4/21,2014
void ordran(int, double[]);// modified 4/21,2014
