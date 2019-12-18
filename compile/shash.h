#ifndef __SHASH_H
#define __SHASH_H

/*
 * Simple string -> fixed-size data hash table.  Note that this data
 * structure accomodates "zero size" data as well, in which case the
 * hash just acts as a table of string entries, and the data pointer
 * returned by shash_index will be the string itself.  shash_set and
 * shash_get should not be used in this case.
 *
 * shash_create   -- create a new shash table
 * shash_pcreate  -- create a new shash table on an existing pool
 * shash_destroy  -- destroy an shash table
 *
 * shash_set      -- set a (key, value) pair in the table
 *                   (overwrites any pair with the same key)
 * shash_get      -- get a value corresponding to the key into a buffer,
 *                   and return it.  If doesn't exist, return null.
 * shash_remove   -- remove an entry corresponding to some key.
 *                   returns 0 on success, -1 if the key didn't exist
 * shash_contains -- check if a key is in the hash table
 * shash_index    -- get a pointer to the value location for a key.
 *                   If no such key, then create a new entry with 0 value.
 *
 * shashc_head    -- Move cursor to the start of the table
 * shashc_next    -- Move the cursor to the next entry
 * shashc_key     -- Look at the key under the cursor
 * shashc_data    -- Look at the data under the cursor
 * shashc_end     -- Return true if the cursor is at end/invalid
 */

#include "mempool.h"

typedef struct shash_struct* shash_t;

typedef struct {
    shash_t hash;
    int bucket_number;
    void* bucket_entry;
} shashc_t;

shash_t  shash_create  (int table_size, int value_size);
shash_t  shash_pcreate (int table_size, int value_size, mempool_t pool);
void     shash_destroy (shash_t hash);
void     shash_set     (shash_t hash, char* key, void* value);
void*    shash_get     (shash_t hash, char* key, void* value);
int      shash_remove  (shash_t hash, char* key);
int      shash_contains(shash_t hash, char* key);
void*    shash_index   (shash_t hash, char* key);

void     shashc_head   (shash_t hash, shashc_t* cursor);
int      shashc_next   (shashc_t* cursor);
char*    shashc_key    (shashc_t* cursor);
void*    shashc_data   (shashc_t* cursor);
int      shashc_end    (shashc_t* cursor);

#endif /* __SHASH_H */
