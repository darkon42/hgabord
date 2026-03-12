
/****************************************************************************/
/*                                                                          */
/*  structpath.h         Paths for mpp files                         */
/*                                                                          */
/****************************************************************************/
/* Previously in help.c */
extern char HELP_DIR[STRING_SIZE]; 

/* Previously in streams_io.c */
extern char MACRO_DIRECTORY[STRING_SIZE];

/*
 * filter path
 */
extern char FltrPath[STRING_SIZE];
/*
 * gabsignal path
 */
extern char SigPath[STRING_SIZE];
/*
 * enviroment table
 */
#ifdef __unix__
    extern char **environ;
#endif


