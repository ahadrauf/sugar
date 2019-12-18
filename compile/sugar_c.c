#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "sugar-lib.h"
#include "lexer.h"
#include "parse.h"
#include "evalexpr.h"
#include "codegen.h"
#include "codeprint.h"
#include "modelmgr.h"

#ifdef USE_XT
#include "xtnetdraw.h"
#endif /* USE_XT */

#ifdef __SUGAR_MEX
#include "codemx.h" 
#else
#include "dcsolve.h"
#include "paramprint.h"
#include "writegeom.h"
#endif


/* If no output file name given, we need to remove a ".net" suffix,
 * if it exists, and tag on a ".m" suffix to get the output file
 * name. 
  */
char* make_outfile_name(char* infile)
{
    static char buf[512];
    int namelen = strlen(infile);

    if (namelen >= sizeof(buf) - 3) {
        return NULL;
    }
    strcpy(buf, infile);
    if (namelen > 4 && strcmp(buf + namelen-4, ".net") == 0)
        strcpy(buf + namelen-4, ".m");
    else
        strcat(buf, ".m");
    return buf;
}

       
#ifndef __SUGAR_MEX

struct compile_flags {
    char* in_file_name;
    char* out_file_name;
    char* uses_file_name;
    char* param_file_name;
    char* geom_file_name;
    shash_t param_bindings;
    int display;
    int do_dc;
    int codegen_flag;
};

int process_compile_flags(int argc, char** argv, struct compile_flags* flags)
{
    flags->in_file_name = NULL;
    flags->out_file_name = NULL;
    flags->uses_file_name = NULL;
    flags->param_file_name = NULL;
    flags->geom_file_name = NULL;
    flags->param_bindings = NULL;
    flags->display = 0;
    flags->do_dc = 0;
    flags->codegen_flag = 0;

    --argc; ++argv;

    if (argc == 0) {
	fprintf(stderr, "sugar-c infile [options]\n");
	fprintf(stderr, "Options are\n");
	fprintf(stderr, "  -o outfile : generate a Matlab script for the netlist\n");
	fprintf(stderr, "  -u uses_file : write out a list of files used in the netlist\n");
	fprintf(stderr, "  -p paramfile : write out a list of parameters used in the netlist\n");
	fprintf(stderr, "  -g geomfile : write a geometry file for display by the Java viewer applet\n");
	fprintf(stderr, "  -d : view the device in an X window\n");
	fprintf(stderr, "  -s : do static analysis\n");
	return -1;
    }
    
    while (argc > 0) {
        char* arg = argv[0];
	++argv; --argc;
	if (*arg == '-') {
	    switch (arg[1]) {
	        case 'o':
		    if (argc < 1) {
			fprintf(stderr, "No output file given");
			return -1;
		    }
		    flags->out_file_name = *argv++;
		    flags->codegen_flag = 1;
		    --argc;
		    break;
		    
		case 'u':
		    if (argc < 1) {
			fprintf(stderr, "No uses output file given");
			return -1;
		    }
		    flags->uses_file_name = *argv++;
		    --argc;
		    break;

		case 'p':
		    if (argc < 1) {
			fprintf(stderr, "No param output file given");
			return -1;
		    }
		    flags->param_file_name = *argv++;
		    --argc;
		    break;

		case 'g':
		    if (argc < 1) {
			fprintf(stderr, "No geom output file given");
			return -1;
		    }
		    flags->geom_file_name = *argv++;
		    flags->codegen_flag = 1;
		    --argc;
		    break;

		case 'd':
		    flags->display = 1;
		    flags->codegen_flag = 1;
		    break;
		    
		case 's':
		    flags->do_dc = 1;
		    flags->codegen_flag = 1;
		    break;
		    
		default:
		    fprintf(stderr, "Unknown flag %s", arg);
		    return -1;
	    }
	} else if (flags->in_file_name == NULL) {
            flags->in_file_name = arg;
	} else {
	    fprintf(stderr, "Cannot compile %s; input file %s already chosen\n",
		    arg, flags->in_file_name);
	    return -1;
	}
    }
    return 0;
}

int main(int argc, char** argv)
{
    extern FILE* yyin;  /* Yacc's input file pointer*/
    extern FILE* yyout; /* Yacc's output file pointer */
    code_t* code;
    netlist_ir* netlist;

    mempool_t parse_pool;
    mempool_t code_pool;

    struct compile_flags flags;
    if (process_compile_flags(argc, argv, &flags) < 0) {
	exit(-1);
    }

    /* Set parser input and output files (plus sanity checks) */

    if (flags.in_file_name == NULL) {
        yyin = stdin;
    } else {
        yyin = fopen(flags.in_file_name, "r");
        if (yyin == NULL) {
            fprintf(stdout, "Cannot open %s\n", flags.in_file_name);
            exit(1);
        }
    }


    /* Parse the input and clean up */

    model_manager_init();
    code_pool = mempool_create(MEMPOOL_DEFAULT_SPAN);
    parse_pool = mempool_create(MEMPOOL_DEFAULT_SPAN);
    code = sugar_parse(parse_pool);
    if (flags.codegen_flag)
        netlist = gen_code(code_pool, code, NULL);
    
    if (no_errors()) {

        if (flags.out_file_name != NULL) {
            yyout = fopen(flags.out_file_name, "w");
            if (yyout == NULL) {
                fprintf(stderr, "Cannot open %s\n", flags.out_file_name);
                exit(1);
            }
	    print_code(yyout, netlist); 
	    fclose(yyout);
	}

	if (flags.uses_file_name != NULL) {
	    FILE* uses_file = fopen(flags.uses_file_name, "w");
	    if (uses_file == NULL) {
		fprintf(stderr, "Cannot open %s\n", flags.uses_file_name);
		exit(1);
            }
	    if (flags.in_file_name) {
		fprintf(uses_file, "%s uses:\n", flags.in_file_name);
            } else {
	        fprintf(uses_file, "uses:\n");
	    }
	    print_files_used(uses_file);
	    fclose(uses_file);
	}
	
	if (flags.param_file_name != NULL) {
	    FILE* param_file = fopen(flags.param_file_name, "w");
	    if (param_file == NULL) {
		fprintf(stderr, "Cannot open %s\n", flags.param_file_name);
		exit(1);
            }
	    print_params(param_file, code);
	    fclose(param_file);
	}
	
	if (flags.geom_file_name != NULL) {
	    FILE* geom_file = fopen(flags.geom_file_name, "w");
	    if (geom_file == NULL) {
		fprintf(stderr, "Cannot open %s\n", flags.geom_file_name);
		exit(1);
            }
	    writegeom(geom_file, netlist, NULL);
	    fclose(geom_file);
	}
	
	/*
        assemble_K(code_pool, netlist);
        assemble_F(code_pool, netlist);
	*/

	if (flags.display) {
#ifdef USE_XT
            display_device(argc, argv, netlist, NULL);
#else
	    fprintf(stderr, "No display capabilities linked in\n");
#endif
	}

	if (flags.do_dc) {
#ifdef USE_XT
	    display_device(argc, argv, netlist, dcsolve(code_pool, netlist));
#else
	    fprintf(stderr, "No display capabilities linked in\n");
#endif
	}
    }

    mempool_destroy(parse_pool);
    mempool_destroy(code_pool);
    model_manager_shutdown();

    if (yyin != stdin)
        fclose(yyin);

    return 0;
}

#else /* __SUGAR_MEX */

static shash_t get_param_bindings(mempool_t parse_pool,
                                  const mxArray* matlab_bindings)
{
    shash_t result = shash_create(127, sizeof(val_t));
    int num_params = mxGetNumberOfFields(matlab_bindings);
    int i;
    
    for (i = 0; i < num_params; ++i) {
        const char* name = mxGetFieldNameByNumber(matlab_bindings, i);
        mxArray* field = mxGetFieldByNumber(matlab_bindings, 0, i);
        val_t value;

        if (mxIsChar(field)) {
            int buflen = mxGetM(field) + 1;
            value.type = 's';
            value.val.s = (char*) mempool_get(parse_pool, buflen);
            mxGetString(field, value.val.s, buflen);
        } else {
            double d = mxGetScalar(field);
            value.type = 'd';
            value.val.d = d;
        }

        shash_set(result, (char*) name, &value);
    }
    return result;
}

void
mexFunction( int nlhs, mxArray *plhs[],
             int nrhs, const mxArray *prhs[] )
{
    char* outfile;
    char* infile;
    FILE* infile_fp;
    extern FILE *yyin;
    extern FILE *yyout;
    netlist_ir* netlist;
    mempool_t parse_pool;
    mempool_t code_pool;
    shash_t params = NULL;

    int namelen;

    if (nrhs == 0) {
        mexErrMsgTxt("Not enough input arguments.\n");
        return;
    }
    if (!mxIsChar(prhs[0]) || (nrhs > 1 && !mxIsChar(prhs[0]))) {
          mexErrMsgTxt("Variable must contain a string.\n");
          return;
    }

    /* Get the name of the input file out of the first argument */
    namelen = mxGetN(prhs[0])+1;
    infile = (char*) mxCalloc(namelen, sizeof(char));
    mxGetString(prhs[0], infile, namelen);

    /* Construct output file name from input name
     */
    outfile = make_outfile_name(infile);
    if (outfile == NULL) {
        mexErrMsgTxt("Name is too long.\n");
        return;
    }

    /* Set parser input and output files (plus sanity checks).
     * Note the use of yyrestart... docs claim that the scanner can become
     * unhappy if you directly set yyin when the scanner is called multiple
     * times
     */
    infile_fp = which_file(infile); 
    if (infile_fp == NULL) {
        char errbuf[256];
        sprintf(errbuf, "Cannot open %s\n", infile);
        mexErrMsgTxt(errbuf);
        return;
    }
    yyrestart(infile_fp);

    if (nlhs < 1) {
        yyout = fopen(outfile, "w");
        if (yyout == NULL) {
            char errbuf[256];
            fclose(infile_fp);
            sprintf(errbuf, "Cannot open %s\n", outfile);
            mexErrMsgTxt(errbuf);
            return;
        }
    }

    /* Parse the input and clean up */

    model_manager_init();
    parse_pool = mempool_create(MEMPOOL_DEFAULT_SPAN);
    code_pool = mempool_create(MEMPOOL_DEFAULT_SPAN);
    if (nrhs > 1) {
        params = get_param_bindings(parse_pool, prhs[1]);
    }
    netlist = gen_code(code_pool, sugar_parse(parse_pool), params);
    mempool_destroy(parse_pool);
    if (no_errors()) {
        if (nlhs < 1)
            print_code(yyout, netlist); 
        else 
            plhs[0] = netlist_to_mx(netlist);
    }
    mempool_destroy(code_pool);
    model_manager_shutdown();

    if (nlhs < 1)
        fclose(yyout);

    fclose(infile_fp);

    if (!no_errors())
        mexErrMsgTxt("There were errors in the netlist.");
}

#endif /* __SUGAR_MEX */
