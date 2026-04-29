/*--------------------------------------------------------------------------*/
/*  MPP gabsignal processing program.                                       */
/* (C) 1993 Copyright New York University, All Rights Reserved.             */
/*                                                                          */
/*  ------------------------------------------------------------            */
/*  Francois Bergeaud, Mike Orszag.                                         */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*  struct_alloc.c   Functions which deal with the memory allocation        */
/*                   of BOOK, WORD and INDEX structures                     */
/*                                                                          */
/*--------------------------------------------------------------------------*/

#include "mpp.h"
#include <stdio.h>
#include <stdlib.h>

/* Fix 9: proper prototypes at file scope replacing all K&R in-body
 * declarations. */
void  init_book(BOOK book);
void  init_word(WORD strt);
void  init_index(INDEX index);
INDEX AllocIndex(void);
INDEX IndexFree(INDEX index);
WORD  WordFree(WORD strt);
WORD  WordListFree(WORD word);

/*--------------------------------------------------------------------------*/
/*
 * AllocBook — allocate and initialise a BOOK structure.
 *
 * Returns the allocated BOOK, or NULL on failure.
 * Caller is responsible for freeing via clear_book() + free().
 *
 * Fix 4: perror() → fprintf+return NULL so the caller receives NULL
 *        instead of crashing on init_book(NULL).
 */
/*--------------------------------------------------------------------------*/
BOOK AllocBook(void)
{
    BOOK book;

    book = (BOOK)malloc(sizeof(struct book));
    /* Fix 4: perror() → fprintf+return NULL */
    if (book == (BOOK)NULL) {
        fprintf(stderr, "AllocBook: malloc failed\n");
        return (BOOK)NULL;
    }

    init_book(book);
    return book;
}

/*--------------------------------------------------------------------------*/
/*
 * init_book — zero-initialise all fields of a BOOK.
 *
 * Fix 4: perror() → fprintf+return so execution stops on NULL input.
 * Fix 11: removed commented-out dead code (book->type = (int)NULL).
 */
/*--------------------------------------------------------------------------*/
void init_book(BOOK book)
{
    /* Fix 4: perror() → fprintf+return */
    if (book == (BOOK)NULL) {
        fprintf(stderr, "init_book(): null argument!\n");
        return;
    }

    book->size     = 0;
    book->energy   = 0.0;
    book->sigen    = 0.0;
    book->type     = 0;       /* Fix 11: removed dead //book->type=(int)NULL */
    book->smin     = 0.0;
    book->smax     = 0.0;
    book->sig_size = 0;
    book->first    = (WORD)NULL;
    book->last     = (WORD)NULL;
}

/*--------------------------------------------------------------------------*/
/*
 * AllocWord — allocate and initialise a WORD structure.
 *
 * Returns the allocated WORD, or NULL on failure.
 * Caller is responsible for freeing via WordFree().
 *
 * Fix 2: if AllocIndex() fails inside init_word, the word struct is now
 *        freed before returning NULL rather than returning a broken WORD.
 * Fix 4: perror() → fprintf+return NULL.
 */
/*--------------------------------------------------------------------------*/
WORD AllocWord(void)
{
    WORD str;

    str = (WORD)malloc(sizeof(struct word));
    /* Fix 4: perror() → fprintf+return NULL */
    if (str == (WORD)NULL) {
        fprintf(stderr, "AllocWord: malloc failed\n");
        return (WORD)NULL;
    }

    str->index = (INDEX)NULL;
    init_word(str);

    /* Fix 2: if init_word failed to allocate the index, free the word
     * struct and return NULL rather than returning a broken WORD. */
    if (str->index == (INDEX)NULL) {
        free(str);
        return (WORD)NULL;
    }

    return str;
}

/*--------------------------------------------------------------------------*/
/*
 * init_word — initialise all fields of a WORD.
 *
 * Fix 5: warning() → fprintf (doesn't stop execution either, but at
 *        least goes to stderr consistently).
 * Fix 10: typo "arugment" corrected to "argument".
 */
/*--------------------------------------------------------------------------*/
void init_word(WORD strt)
{
    if (strt == (WORD)NULL) {
        /* Fix 5/10: warning() + typo → fprintf */
        fprintf(stderr, "init_word(): null argument!\n");
        return;
    }

    strt->coeff  = 0.0;
    strt->value  = 0.0;
    strt->coeff1 = 0.0;
    strt->value1 = 0.0;
    strt->next   = (WORD)NULL;
    strt->status = KEPT;

    if (strt->index == (INDEX)NULL)
        strt->index = AllocIndex();
    else
        init_index(strt->index);
}

/*--------------------------------------------------------------------------*/
/*
 * AllocIndex — allocate and initialise an INDEX structure.
 *
 * Returns the allocated INDEX, or NULL on failure.
 * Caller is responsible for freeing via IndexFree().
 *
 * Fix 4: perror() → fprintf+return NULL.
 */
/*--------------------------------------------------------------------------*/
INDEX AllocIndex(void)
{
    INDEX indx;

    indx = (INDEX)malloc(sizeof(struct index));
    /* Fix 4: perror() → fprintf+return NULL */
    if (indx == (INDEX)NULL) {
        fprintf(stderr, "AllocIndex: malloc failed\n");
        return (INDEX)NULL;
    }

    init_index(indx);
    return indx;
}

/*--------------------------------------------------------------------------*/
/*
 * init_index — zero-initialise all fields of an INDEX.
 *
 * Fix 5: warning() → fprintf+return.
 */
/*--------------------------------------------------------------------------*/
void init_index(INDEX index)
{
    if (index == (INDEX)NULL) {
        /* Fix 5: warning() → fprintf+return */
        fprintf(stderr, "init_index(): null argument!\n");
        return;
    }

    index->id       = 0.0;
    index->octave   = 0.0;
    index->position = 0.0;
    index->phase    = 0.0;
}

/*--------------------------------------------------------------------------*/
/*
 * clear_book — free all WORDs in the book and re-initialise its fields.
 *
 * Fix 5/7: warning() → fprintf+return so execution stops on NULL input
 *          instead of crashing on book->first below.
 */
/*--------------------------------------------------------------------------*/
void clear_book(BOOK book)
{
    if (book == (BOOK)NULL) {
        /* Fix 5/7: warning() → fprintf+return */
        fprintf(stderr, "clear_book(): null argument!\n");
        return;
    }

    WordListFree(book->first);
    init_book(book);
}

/*--------------------------------------------------------------------------*/
/*
 * WordFree — free a single WORD and its INDEX.
 *
 * Returns NULL so the caller can write: word = WordFree(word).
 *
 * Fix 5/6: warning() → fprintf+return NULL so execution stops on NULL
 *          input instead of crashing on strt->index below.
 */
/*--------------------------------------------------------------------------*/
WORD WordFree(WORD strt)
{
    if (strt == (WORD)NULL) {
        /* Fix 5/6: warning() → fprintf+return NULL */
        fprintf(stderr, "WordFree(): null argument!\n");
        return (WORD)NULL;
    }

    IndexFree(strt->index);
    free(strt);
    return (WORD)NULL;
}

/*--------------------------------------------------------------------------*/
/*
 * WordListFree — free an entire linked list of WORDs.
 *
 * Returns NULL so the caller can write: book->first = WordListFree(book->first).
 *
 * Fix 3: was recursive — O(n) stack frames, stack overflow risk for large
 *        books. Replaced with an iterative loop using O(1) stack.
 *        The next pointer is saved before WordFree() is called because
 *        WordFree() frees the node, after which reading strt->next would
 *        be undefined behaviour.
 */
/*--------------------------------------------------------------------------*/
WORD WordListFree(WORD word)
{
    WORD next;

    while (word != (WORD)NULL) {
        next = word->next;   /* save before freeing */
        WordFree(word);
        word = next;
    }
    return (WORD)NULL;
}

/*--------------------------------------------------------------------------*/
/*
 * IndexFree — free an INDEX structure.
 *
 * Returns NULL so the caller can write: index = IndexFree(index).
 *
 * Fix 5/8: warning() → fprintf+return NULL. free(NULL) is safe per the
 *          C standard but the warning was misleading — NULL is now just
 *          a no-op with a logged message.
 */
/*--------------------------------------------------------------------------*/
INDEX IndexFree(INDEX index)
{
    if (index == (INDEX)NULL) {
        /* Fix 5/8: warning() → fprintf+return NULL */
        fprintf(stderr, "IndexFree(): null argument!\n");
        return (INDEX)NULL;
    }

    free(index);
    return (INDEX)NULL;
}

/*--------------------------------------------------------------------------*/
/*
 * IndexAlloc — allocate an INDEX without initialisation.
 *
 * Returns the allocated INDEX, or NULL on failure.
 * Prefer AllocIndex() which also initialises the fields.
 *
 * Fix 4: perror() → fprintf+return NULL.
 */
/*--------------------------------------------------------------------------*/
INDEX IndexAlloc(void)
{
    INDEX index;

    index = (INDEX)malloc(sizeof(struct index));
    /* Fix 4: perror() → fprintf+return NULL */
    if (index == (INDEX)NULL) {
        fprintf(stderr, "IndexAlloc: malloc failed\n");
        return (INDEX)NULL;
    }

    return index;
}

/*--------------------------------------------------------------------------*/
/*
 * end of struct_alloc.c
 */
/*--------------------------------------------------------------------------*/