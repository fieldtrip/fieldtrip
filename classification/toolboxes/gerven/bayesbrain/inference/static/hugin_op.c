/*
 * creates a junction tree using the HUGIN library.
 *
 * compile with:
 * mex -I/usr/local/hugin/include -L/usr/local/hugin/lib -lhugin -lm -lz -lcrypt -lpthread hugin_op.c
 * path should be defined in LD_LIBRARY_PATH environment variable
 *
 * implementation of switchyard to make dom persistent
 *
 */

# include "mex.h"
# include "hugin.h"

# define max(a,b) (a>b ? a : b)

/* persistent variable */
h_domain_t dom;    
    

/* initialize junction tree */
void init(const mxArray *prhs)
{
    char *input_buf;    
    mwSize buflen;

    /* input must be a string */
    if ( mxIsChar(prhs) != 1)
        mexErrMsgTxt("Input must be a string.");
    
        /* get the length of the input string */
    buflen = (mxGetM(prhs) * mxGetN(prhs)) + 1;
    
        /* copy the string data from prhs[2] into a C string input_ buf.    */
    input_buf = mxArrayToString(prhs);
    
    if(input_buf == NULL)
        mexErrMsgTxt("Could not convert input to string.");
    
    dom = h_net_parse_domain(input_buf, NULL,NULL);
    
    h_domain_triangulate(dom, h_tm_clique_weight);
    
    h_domain_compile(dom);

    mxFree(input_buf);

}

/* enter evidence */
void enter_evidence(const mxArray *pevid)
{
   
    h_node_t nd;
    size_t state;
    
    double *xValues;
    
    int i;
    int nval;
    
    char str[10];
        
    xValues = mxGetPr(pevid);
    
    nval = max(mxGetN(pevid),mxGetM(pevid));
    
    for(i=0;i<nval;i++)
    {
        if (!mxIsNaN(xValues[i]))
        {
            sprintf(str,"C%d",i+1);
            
            nd = h_domain_get_node_by_name(dom, str);
            
            if (h_node_get_kind(nd) == h_kind_continuous)
                h_node_enter_value(nd, (h_double_t)(xValues[i])); /* continuous */
            else
                h_node_select_state(nd, (size_t)(xValues[i]-1)); /* discrete */
            
        }
    }
         
    h_domain_propagate(dom, h_equilibrium_sum, h_mode_normal);
  
}


/* marginalize */
mxArray * marginalize(const mxArray *prhs)
{
       
    int i,K;
    
    mxArray *plhs;
    double *outArray;
    
    double *xValues;
    char str[10];
    
    xValues = mxGetPr(prhs);
    
    sprintf(str,"C%d",(int)(xValues[0]));
        
    h_node_t nd = h_domain_get_node_by_name(dom, str);
        
    
    
    if (h_node_get_kind(nd) == h_kind_continuous)
    {
        plhs = mxCreateDoubleMatrix(1, 2, mxREAL);
        
        /* Get a pointer to the data space in our newly allocated memory */
        outArray = mxGetPr(plhs);

        outArray[0] = h_node_get_mean(nd);
        
        outArray[1] = h_node_get_variance(nd);
    }
    else
    {
        K = h_node_get_number_of_states(nd);
    
        plhs = mxCreateDoubleMatrix(1, K, mxREAL);
        
        /* Get a pointer to the data space in our newly allocated memory */
        outArray = mxGetPr(plhs);
        
        /* return discrete data */
        for(i=0;i<K;i++)
        {
            outArray[i] = h_node_get_belief(nd,i);
        }
    }
    
    return plhs;
}

/* switchyard */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char *option; 
    mwSize buflen;
    
    /* input must be a string */
    if ( mxIsChar(prhs[0]) != 1)
      mexErrMsgTxt("Input must be a string.");
    
    /* get the length of the input string */
    buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;

    /* copy the string data from prhs[0] into a C string input_ buf.    */
    option = mxArrayToString(prhs[0]);

    if(option == NULL)
      mexErrMsgTxt("Could not convert input to string.");
    
    if (strcmp("init",option) == 0)
        
        init(prhs[1]);
  
    else if (strcmp("enter",option) == 0)
  
        enter_evidence(prhs[1]);
    
    else if (strcmp("marg",option) == 0)
        
        plhs[0] = marginalize(prhs[1]);
    
    else 
        mexPrintf("unknown option\n");
    
    mxFree(option);
}
