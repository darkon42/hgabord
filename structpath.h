/****************************************************************************/
/*                                                                          */
/*  structpath.h    Paths for mpp files                                     */
/*                                                                          */
/****************************************************************************/

/* Fix 14: include guard added — was missing. */
#ifndef STRUCTPATH_H
#define STRUCTPATH_H

/* Fix 15: STRING_SIZE is defined in mpp.h / gabsignals.h. Guard it here
 * so this header is safe to include even if mpp.h has not been included
 * first, and to avoid redefinition warnings if it has. */
#ifndef STRING_SIZE
#define STRING_SIZE 200
#endif

/* Previously in help.c */
extern char HELP_DIR[STRING_SIZE];

/* Previously in streams_io.c */
extern char MACRO_DIRECTORY[STRING_SIZE];

/*
 * filter path
 */
extern char FltrPath[STRING_SIZE];

/*
 * signal path
 */
extern char SigPath[STRING_SIZE];

/*
 * environment table
 */
#ifdef __unix__
    extern char **environ;
#endif

/*
 * end of structpath.h
 */

#endif /* STRUCTPATH_H */