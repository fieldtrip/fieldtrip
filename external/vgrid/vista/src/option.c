/*
 *  Copyright 1994 University of British Columbia
 *
 *  Permission to use, copy, modify, distribute, and sell this software and its
 *  documentation for any purpose is hereby granted without fee, provided that
 *  the above copyright notice appears in all copies and that both that
 *  copyright notice and this permission notice appear in supporting
 *  documentation. UBC makes no representations about the suitability of this
 *  software for any purpose. It is provided "as is" without express or
 *  implied warranty.
 *
 *  Author: Arthur Pope, UBC Laboratory for Computational Intelligence
 */

/*
 *  Remodeled for SimBio by F. Kruggel (kruggel@cns.mpg.de) - 15/07/00
 */

#include <vista.h>
#define __USE_BSD		/* to define S_ISSOCK & S_ISLNK */
#include <sys/stat.h>

static int ParseArgValues (int *arg, int argc,
			   char **argv, VOptionDescRec *opt)
{
  VArgVector *vec = (VArgVector *) opt->value;
  VPointer values;
  int nvalues = 0;
  VDictEntry *dict = opt->dict;
  char *cp;

  /* Locate the place we're to store the argument values: */
  if (opt->number == 0) {
    /* If a variable number of arguments is expected, allocate storage
       to hold them: */
    vec->vector = values = VMalloc ((argc - *arg) * VRepnSize (opt->repn));
  } else values = opt->value;

  /* If no dictionary is specified for a boolean-valued option, use the
     default one of true, false, yes, no... */
  if (opt->repn == VBooleanRepn && ! dict)
    dict = VBooleanDict;

  /* Parse argument values until we've reached the required number,
     we've run out of entered arguments, or we encounter one that
     is ill-formed: */
  while ((opt->number == 0 || nvalues < opt->number) && (*arg < argc)) {
    cp = argv[*arg];

    /* Special treatment for string-valued options: */
    if (opt->repn == VStringRepn) {

      /* An argument of the form -string is not interpreted as a
	 string value: */
      if (cp[0] == '-' && cp[1] != 0 && cp[1] != '-')
	break;

      /* An argument of the form --string is interpreted as string: */
      if (cp[0] == '-' && cp[1] == '-' && cp[2] != 0)
	cp += 2;
    }

    /* Convert the argument to the specified internal form: */
    if (! VDecodeAttrValue (cp, dict, opt->repn, values))
      break;
    nvalues++;
    values = (VPointer) ((char *) values + VRepnSize (opt->repn));
    argv[(*arg)++] = NULL;
  }

  /* Special treatment of boolean-valued options: if the option has just
     one value associated with it then treat -option <other options> like
     -option true <other options>: */
  if (opt->repn == VBooleanRepn && opt->number == 1 && nvalues == 0) {
    * (VBoolean *) opt->value = TRUE;
    nvalues = 1;
  }

  return nvalues;
}

VBoolean VIdentifyFiles (int noptions, VOptionDescRec options[],
			 VStringConst keyword, int *argc, char **argv, int fd)
{
  VOptionDescRec *opt;
  VArgVector *vec;
  VBoolean switch_found = FALSE;
  VStringConst *values;
  struct stat stat_buf;
  int i, j, n;

  /* A -in or -out switch on the command line may specify one or more
     files. Look for the switch's entry in the option table supplied: */
  for (i = 0, opt = options; i < noptions; i++, opt++)
    if (strcmp (keyword, opt->keyword) == 0)
      break;
  if (i == noptions)
    VError ("VIdentifyFiles: Option -%s not defined in option table",
	    keyword);

  /* The option table entry must specify a location for storing filenames: */
  if (! opt->value)
    VError ("VIdentifyFiles: No value storage for option -%s", keyword);
  vec = ((VArgVector *) opt->value);

  /* Whether or not the option was present on the command line should have
     been recorded in a "found" flag: */
  if (opt->found && opt->found != VRequiredOpt && opt->found != VOptionalOpt)
    switch_found = *(opt->found);
  else VError ("VIdentifyFiles: No \"found\" flag for option -%s", keyword);
	
  /* If a -in or -out switch was specified, it supplies the files' names: */
  if (switch_found)
    return TRUE;

  /* If no switch was specified but some extra arguments were supplied
     with the command, those arguments should be the files' names: */
  if (*argc > 1) {

    /* If any number are permitted, allocate storage for all of them: */
    if (opt->number == 0) {
      vec->vector = VMalloc ((*argc - 1) * sizeof (VStringConst));
      values = (VStringConst *) vec->vector;
    } else values = (VStringConst *) opt->value;

    /* Accept arguments as filenames provided they don't begin with
       a single -: */
    for (i = j = 1, n = 0; i < *argc; i++)
      if ((argv[i][0] == '-' && argv[i][1] != '-') ||
	  (opt->number && n == opt->number))
	argv[j++] = argv[i];		/* not a filename */
      else if (argv[i][0] == '-' && argv[i][1] == '-')
	values[n++] = argv[i] + 1;	/* filename with -- prefix */
      else values[n++] = argv[i];		/* filename */

    /* Ensure that a sufficient number were supplied: */
    if (n) {
      if (opt->number) {
	if (n != opt->number) {
	  if (opt->number == 1)	goto WrongNumber1;
	  else goto WrongNumberN;
        };
      } else vec->number = n;
      *argc = j;
      return TRUE;
    }
  }

  /* That leaves only an open file descriptor as the file. In that case
     there can be at most one such file: */
  if (opt->number > 1) {
  WrongNumberN:
    fprintf (stderr, "%s: %d files must be specified by -%s or "
	     "extra command arguments.\n", argv[0], opt->number, keyword);
    return FALSE;
  }

  /* If a particular file descriptor is specified by fd, ensure that it
     is connected to a pipe or file (not the terminal): */
  if (fd >= 0) {
    if (fstat (fd, & stat_buf) == -1)
      VSystemError ("Failed to fstat() fd %d", fd);
    if (! S_ISREG (stat_buf.st_mode) && ! S_ISSOCK (stat_buf.st_mode) &&
	! S_ISLNK (stat_buf.st_mode) && ! S_ISFIFO (stat_buf.st_mode)) {
	WrongNumber1:
      fprintf (stderr, "%s: No file specified by -%s, extra command "
	       "argument, or %c.\n", argv[0], keyword, "<>"[fd]);
      return FALSE;
    }
  }

  /* Record the file as being "-": */
  if (opt->number == 0) {
    vec->number = 1;
    vec->vector = VMalloc (sizeof (VStringConst));
    ((VStringConst *) vec->vector)[0] = "-";
  } else * ((VStringConst *) opt->value) = "-";
  return TRUE;
}

VBoolean VParseCommand (int noptions, VOptionDescRec options[],
			int *argc, char **argv)
{
  int arg, nvalues, i, j;
  char *cp;
  VBoolean *opts_seen, result = TRUE;
  VOptionDescRec *opt, *opt_t;

  /* Allocate storage for a set of flags indicating which
     arguments have been seen: */
  opts_seen = VCalloc (noptions, sizeof (VBoolean));

  /* Initialize any "found" flags to false, and the number field of any
     VArgVector values to zero: */
  for (opt = options + noptions - 1; opt >= options; opt--) {
    if (opt->found)
      *opt->found = FALSE;
    if (opt->number == 0 && opt->value)
      ((VArgVector *) opt->value)->number = 0;
  }

  /* For each argument supplied with the command: */
  for (arg = 1; arg < *argc; ) {
    cp = argv[arg++];

    /* If it doesn't start with - it can't be an option: */
    if (cp[0] != '-' || cp[1] == 0)
      continue;

    /* Check for -help: */
    if (strcmp (cp + 1, "help") == 0) {

      /* If found, return FALSE to force printing of usage info: */
      *argc = 1;
      return FALSE;
    }

    /* Look up the argument in the list of options: */
    i = strlen (cp + 1);
    opt = NULL;
    for (opt_t = options + noptions - 1; opt_t >= options; opt_t--) {
      if (strncmp (cp + 1, opt_t->keyword, i) != 0)
	continue;	/* not this one */
      if (i == strlen (opt_t->keyword)) {
	opt = opt_t;
	break;		/* an exact match */
      }
      if (opt)
	goto NextArg;	/* already matched another prefix */
      opt = opt_t;	/* note a prefix match */
    }

    /* If the argument isn't recognized, skip it: */
    if (! opt)
      goto NextArg;	/* not recognized */

    /* Remove it from the list of command arguments: */
    argv[arg - 1] = 0;

    /* Ensure that the option has not already been seen: */
    if (opts_seen[opt - options]) {
      fprintf (stderr,
	       "%s: Duplicate -%s option; ignoring all but last.\n",
	       argv[0], opt->keyword);

      /* If it has been seen, delete its previous value: */
      if (opt->number == 0) {
	VFree (((VArgVector *) opt->value)->vector);
	((VArgVector *) opt->value)->number = 0;
      }
    } else opts_seen[opt - options] = TRUE;

    /* Swallow any value(s) that follow: */
    switch (opt->repn) {

    case VBitRepn:
    case VUByteRepn:
    case VSByteRepn:
    case VShortRepn:
    case VLongRepn:
    case VFloatRepn:
    case VDoubleRepn:
    case VBooleanRepn:
    case VStringRepn:
      nvalues = ParseArgValues (& arg, *argc, argv, opt);
      break;

    default:
      VError ("Parsing of command options with %s values "
	      "is not implemented", VRepnName (opt->repn));
      nvalues = 0;	/* to quiet lint */
    }

    /* Ensure that the expected number of arguments was found: */
    if (opt->number && nvalues != opt->number) {

      /* Either we encountered an argument we couldn't parse, or
	 we used up all arguments before finding the expected number
	 of them: */
      fprintf (stderr, "%s: Option -%s ", argv[0], opt->keyword);
      if (arg < *argc)
	fprintf (stderr, "has incorrect value %s.\n", argv[arg]);
      else if (opt->number > 1)
	fprintf (stderr, "requires %d values; found only %d.\n",
		 opt->number, nvalues);
      else fprintf (stderr, "requires a value.\n");
      result = FALSE;
      break;
    }
    if (opt->number == 0)
      ((VArgVector *) opt->value)->number = nvalues;

    /* Note that a value was successfully obtained for this option: */
    if (opt->found)
      *(opt->found) = TRUE;

  NextArg: ;
  }

  /* Ensure that each mandatory option was seen: */
  for (i = 0; i < noptions; i++)
    if (options[i].found == VRequiredOpt && ! opts_seen[i]) {
      fprintf (stderr, "%s: Option -%s must be specified.\n",
	       argv[0], options[i].keyword);
      result = FALSE;
    }
  VFree ((VPointer) opts_seen);

  /* Squeeze together the remaining arguments in argv: */
  for (i = j = 1; i < *argc; i++)
    if (argv[i])
      argv[j++] = argv[i];
  *argc = j;

  return result;
}

void VParseFilterCmd (int noptions, VOptionDescRec opts[],
		      int argc, char **argv, FILE **inp, FILE **outp)
{
  VStringConst in_file = 0, out_file = 0;
  VBoolean in_found = 0, out_found = 0;
  VOptionDescRec io_opts[] = {
    { "in", VStringRepn, 1, & in_file, & in_found, NULL, "Input file" },
    { "out", VStringRepn, 1, & out_file, & out_found, NULL, "Output file" }
  };
  int i, n;
  VOptionDescRec options[100];
    
  /* Check that number of options will not overflow the options array. */
  if (noptions >= 98) {
    VWarning("VParseFilterCmd: Too many options allowed");
    noptions = 98;
  }

  /* Copy options into a new list and add the "in" and "out" options. */
  n = 0;
  if (inp)
    options[n++] = io_opts[0];
  if (outp)
    options[n++] = io_opts[1];
  for (i = 0; i < noptions; i++, n++)
    options[n] = opts[i];
    
  /* Parse command line arguments and identify the input and output files: */
  if (! VParseCommand (n, options, & argc, argv) ||
      (inp && ! VIdentifyFiles (n, options, "in", & argc, argv, 0)) ||
      (outp && ! VIdentifyFiles (n, options, "out", & argc, argv, 1)))
    goto Usage;

  /* Any remaining unparsed arguments are erroneous: */
  if (argc > 1) {
    VReportBadArgs (argc, argv);
  Usage:  VReportUsage (argv[0], n, options,
			inp ? (outp ? "[infile] [outfile]" : "[infile]") :
			(outp ? "[outfile]" : NULL));
    exit (EXIT_FAILURE);
  }

  /* Open the input and output files: */
  if (inp)
    *inp = VOpenInputFile (in_file, TRUE);
  if (outp)
    *outp = VOpenOutputFile (out_file, TRUE);
}
    
void VPrintOptions (FILE *f, int noptions, VOptionDescRec options[])
{
  int i;

#define commentColumn 40

  for (i = 0; i < noptions; i++, options++) {
    fputs ("  ", f);
    if (options->blurb)
      fprintf (f, "%s ", options->blurb);
    fprintf (f, "(%s): ", options->keyword);
    VPrintOptionValue (f, options);
    fputc ('\n', f);
  }

#undef commentColumn
}

int VPrintOptionValue (FILE *f, VOptionDescRec *option)
{
  int n, i, col = 0;
  char *vp;
  VDictEntry *dict;
  VLong ivalue;
  VDouble fvalue = 0.0;
  VStringConst svalue;

  if (option->number == 0) {
    n = ((VArgVector *) option->value)->number;
    vp = (char *) ((VArgVector *) option->value)->vector;
  } else {
    n = option->number;
    vp = (char *) option->value;
  }
  for (i = 0; i < n; i++, vp += VRepnSize (option->repn)) {
    if (i > 0)
      fputc (' ', f);
    switch (option->repn) {
	    
    case VBitRepn:
      ivalue = * (VBit *) vp;
      goto PrintLong;
	    
    case VUByteRepn:
      ivalue = * (VUByte *) vp;
      goto PrintLong;
	    
    case VSByteRepn:
      ivalue = * (VSByte *) vp;
      goto PrintLong;
	    
    case VShortRepn:
      ivalue = * (VShort *) vp;
      goto PrintLong;
	    
    case VLongRepn:
      ivalue = * (VLong *) vp;
    PrintLong:  if (option->dict &&
		    (dict = VLookupDictValue (option->dict, VLongRepn, ivalue)))
      col += fprintf (f, "%s", dict->keyword);
    else col += fprintf (f, "%d", ivalue);
      break;
	    
    case VFloatRepn:
      fvalue = * (VFloat *) vp;
      goto PrintDbl;
	    
    case VDoubleRepn:
      fvalue = * (VDouble *) vp;
    PrintDbl:   if (option->dict &&
		    (dict = VLookupDictValue (option->dict, VDoubleRepn, fvalue)))
      col += fprintf (f, "%s", dict->keyword);
    else col += fprintf (f, "%g", fvalue);
      break;
	    
    case VBooleanRepn:
      col += fprintf (f, "%s", * (VBoolean *) vp ? "true" : "false");
      break;
	    
    case VStringRepn:
      svalue = * (VString *) vp;
      if (! svalue)
	svalue = "(none)";
      else if (option->dict &&
	       (dict = VLookupDictValue (option->dict, VStringRepn, svalue)))
	svalue = dict->keyword;
      col += fprintf (f, "%s", svalue);
      break;

    default:
      break;
    }
  }

  return col;
}

void VReportBadArgs (int argc, char **argv)
{
  int i;

  if (argc > 1) {
    fprintf (stderr, "%s: Unrecognized argument%s:", argv[0],
	     (argc == 2) ? "" : "s");
    for (i = 1; i < argc; i++)
      fprintf (stderr, " %s", argv[i]);
    fputc ('\n', stderr);
  }
}

void VReportUsage (VStringConst program, int noptions,
		   VOptionDescRec options[], VStringConst other_args)
{
  fprintf (stderr, "\nUsage: %s <options>", program);
  if (other_args)
    fprintf (stderr, " %s", other_args);
  fprintf (stderr, ", where <options> includes:\n");
  fprintf (stderr, "    -help\n\tPrints usage information.\n");

  /* Print a line describing each possible option: */
  VReportValidOptions (noptions, options);
}

void VReportValidOptions (int noptions, VOptionDescRec options[])
{
  int i, j;
  VDictEntry *dict;
  VStringConst cp;

  /* Print a line describing each possible option: */
  for (i = 0; i < noptions; i++, options++) {

    /* Print the name of the option: */
    fprintf (stderr, "    -%s ", options->keyword);

    /* Get a name for the type of value needed: */
    if (VIsIntegerRepn (options->repn))
      cp = "integer";
    else if (VIsFloatPtRepn (options->repn))
      cp = "number";
    else if (options->repn == VBooleanRepn)
      cp = "boolean";
    else cp = "string";

    /* Print possible option values: */
    if (options->dict && options->number == 1) {
      for (dict = options->dict; dict->keyword; dict++) {
	if (dict > options->dict)
	  fputs (" | ", stderr);
	fputs (dict->keyword, stderr);
      }
    } else if (options->repn == VBooleanRepn && options->number == 1) {
      fputs ("[ true | false ]", stderr);
    } else if (options->number <= 1) {
      fprintf (stderr, "<%s>", options->dict ? "keyword" : cp);
      if (options->number == 0)
	fputs (" ...", stderr);
    } else {
      for (j = 1; j <= options->number; j++) {
	if (j > 1)
	  fputc (' ', stderr);
	fprintf (stderr, "<%s%d>", options->dict ? "keyword" : cp, j);
      }
    }
    if (options->dict && options->number != 1) {
      fputs ("\n\t<keyword> is ", stderr);
      for (dict = options->dict; dict->keyword; dict++) {
	if (dict > options->dict)
	  fputs (" | ", stderr);
	fputs (dict->keyword, stderr);
      }
    }
    fputs ("\n\t", stderr);

    /* Print any blurb available, and default value(s): */
    if (options->blurb)
      fprintf (stderr, "%s. ", options->blurb);
    if (options->found == VRequiredOpt)
      fputs ("Required.", stderr);
    else if (options->number > 0) {
      fprintf (stderr, "Default: ");
      VPrintOptionValue (stderr, options);
    }
    fputc ('\n', stderr);
  }
}
