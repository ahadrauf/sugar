#include <assert.h>
#include <string.h>

#include "sugar-lib.h"  /* Used only to override standard mem ops, printf */
#include "mempool.h"
#include "shash.h"

/* typedef struct {
 *     shash_t hash;
 *     int bucket_number;
 *     void* bucket_entry;
 * } shashc_t;
 *
 * For the "end of hash" case, bucket_number = -1 and bucket_entry == NULL.
 */

/* struct shash_bucket {
 *     void* next;
 *     char value[value_size]; // correctly aligned
 *     char key[...];          // not necessarily aligned
 * };
 */

/* Accessor macros for hash buckets */
#define bucket_next(bucketp)  ( *((void**) (bucketp)) )
#define bucket_value(bucketp) (  ((char*) (bucketp)) + align_size )
#define bucket_key(bucketp)   (  bucket_value(bucketp) + hash->value_size )
#define bucket_size(key)      ( align_size + hash->value_size + strlen(key)+1 )

#define SHASH_MAGIC_NUMBER 0x2673
#define ASSERT_VALID_HASH assert(hash != NULL && \
		                 hash->magic_number == SHASH_MAGIC_NUMBER)

struct shash_struct {
    int magic_number;
    mempool_t pool;     /* Pool to use for allocations */
    int table_size;     /* Number of buckets in table  */
    int value_size;     /* Size of data value (may == 0 for string table) */
    int own_pool;       /* True if pool allocated during hash construction */
    void** table;       /* Array of buckets */
};

/* Compute the hash function for key (modulo table_size)
 */
static int compute_hash(char* key, int table_size)
{
    int hash_index = 0;
    while (*key != '\0') {
        hash_index = ((hash_index << 8) + (int) *key) % table_size;
        ++key;
    }
    return hash_index;
}

/* Create a new hash table with its own pool for internal allocations
 */
shash_t shash_create(int table_size, int value_size)
{
    mempool_t pool = mempool_create(MEMPOOL_DEFAULT_SPAN);
    shash_t hash = shash_pcreate(table_size, value_size, pool);
    hash->own_pool = 1;
    return hash;
}

/* Create a new hash table on an existing pool
 */
shash_t shash_pcreate(int table_size, int value_size, mempool_t pool)
{
    shash_t hash;

    assert(table_size > 0 && value_size >= 0);

    hash = (shash_t) mempool_get(pool, sizeof(struct shash_struct));

    hash->magic_number = SHASH_MAGIC_NUMBER;
    hash->pool       = pool;
    hash->table_size = table_size;
    hash->value_size = value_size;
    hash->table      = mempool_get(pool, table_size * sizeof(void*));
    hash->own_pool   = 0;
    memset(hash->table, 0, table_size * sizeof(void*));
    return hash;
}

/* Destroy a hash table
 */
void shash_destroy(shash_t hash)
{
    ASSERT_VALID_HASH;
    hash->magic_number = 0;
    if (hash->own_pool)
        mempool_destroy(hash->pool);
}

/* Copy a value into the entry for key.  Overwrite any existing entry. 
 */
void shash_set(shash_t hash, char* key, void* value)
{
    int hash_index;
    void* hash_bucket;

    ASSERT_VALID_HASH;
    assert(key != NULL && (hash->value_size == 0 || value != NULL));

    hash_index = compute_hash(key, hash->table_size);
    hash_bucket = hash->table[hash_index];

    /* Replace the current value if the key is already there */
    while (hash_bucket != NULL) {
        if (strcmp(bucket_key(hash_bucket), key) == 0) {
            if (hash->value_size > 0)
                memcpy(bucket_value(hash_bucket), value, hash->value_size);
            return;
        }
        hash_bucket = bucket_next(hash_bucket);
    }

    /* Otherwise, add a new bucket to the head of the list */
    hash_bucket = mempool_get(hash->pool, bucket_size(key));
    bucket_next(hash_bucket)  = hash->table[hash_index];
    if (hash->value_size > 0)
       memcpy(bucket_value(hash_bucket), value, hash->value_size);
    strcpy(bucket_key(hash_bucket), key);
    hash->table[hash_index] = hash_bucket;
}

/* Save the entry for key into the passed buffer.
 * If no such entry, return NULL; else return pointer to value.
 */
void* shash_get(shash_t hash, char* key, void* value)
{
    int hash_index;
    void* hash_bucket;

    ASSERT_VALID_HASH;
    assert(key != NULL && (hash->value_size == 0 || value != NULL));

    hash_index = compute_hash(key, hash->table_size);
    hash_bucket = hash->table[hash_index];

    /* Return the key if it is there */
    while (hash_bucket != NULL) {
        if (strcmp(bucket_key(hash_bucket), key) == 0) {
            if (hash->value_size > 0)
                memcpy(value, bucket_value(hash_bucket), hash->value_size);
            return value;
        }
        hash_bucket = bucket_next(hash_bucket);
    }

    return NULL;
}

/* Remove the entry for key from the table.  
 * Return 0 on success, < 0 on failure.
 */
int shash_remove(shash_t hash, char* key)
{
    int hash_index;
    void** hash_bucketp;

    ASSERT_VALID_HASH;
    assert(key != NULL);

    hash_index = compute_hash(key, hash->table_size);
    hash_bucketp = &(hash->table[hash_index]);

    /* Remove the key if it is there */
    while (*hash_bucketp != NULL) {
        if (strcmp(bucket_key(*hash_bucketp), key) == 0) {
            *hash_bucketp = bucket_next(*hash_bucketp);
            return 0;
        }
        hash_bucketp = (void**) *hash_bucketp;
    }

    return -1;
}

/* Return true iff hash contains the given key
 */
int shash_contains(shash_t hash, char* key)
{
    int hash_index;
    void* hash_bucket;

    ASSERT_VALID_HASH;
    assert(key != NULL);

    hash_index = compute_hash(key, hash->table_size);
    hash_bucket = hash->table[hash_index];

    /* Return the key if it is there */
    while (hash_bucket != NULL) {
        if (strcmp(bucket_key(hash_bucket), key) == 0) {
            return 1;
        }
        hash_bucket = bucket_next(hash_bucket);
    }

    return 0;
}

/* Return a pointer to the value slot associated with the bucket for key.
 * If no such entry exists, create one.
 */
void* shash_index(shash_t hash, char* key)
{
    int hash_index;
    void* hash_bucket;

    ASSERT_VALID_HASH;
    assert(key != NULL);

    hash_index = compute_hash(key, hash->table_size);
    hash_bucket = hash->table[hash_index];

    /* Return the key if it is there */
    while (hash_bucket != NULL) {
        if (strcmp(bucket_key(hash_bucket), key) == 0) {
            return bucket_value(hash_bucket);
        }
        hash_bucket = bucket_next(hash_bucket);
    }

    /* Otherwise, create and zero out a new element */
    hash_bucket = mempool_get(hash->pool, bucket_size(key));
    bucket_next(hash_bucket) = hash->table[hash_index];
    memset(bucket_value(hash_bucket), 0, hash->value_size);
    strcpy(bucket_key(hash_bucket), key);
    hash->table[hash_index] = hash_bucket;

    return bucket_value(hash_bucket);
}

/* Set a cursor to the first entry in the hash table.  If the
 * table is entry, the cursor is set to the "end" state.
 */
void shashc_head(shash_t hash, shashc_t* cursor)
{
    int bucket_index;

    ASSERT_VALID_HASH;
    assert(cursor != NULL);

    cursor->hash = hash;

    /* Find first entry, if it exists */
    for (bucket_index = 0; bucket_index < hash->table_size; ++bucket_index) {
        if (hash->table[bucket_index] != NULL) {
            cursor->bucket_number = bucket_index;
            cursor->bucket_entry = hash->table[bucket_index];
            return;
        }
    }

    /* Otherwise, set to end/invalid entry */ 
    cursor->bucket_number = -1;
    cursor->bucket_entry = NULL;
}

/* Common sanity check code for the cursor routines:
 * ensure that the cursor is valid and points to a valid hash,
 * and set the variable "hash" to cursor->hash.
 */
#define SET_VALID_HASH \
    shash_t hash; \
    assert(cursor != NULL); \
    hash = cursor->hash; \
    ASSERT_VALID_HASH

/* Move the cursor to the next entry in the table.  Return true iff
 * cursor has moved to a non-end entry.
 */
int shashc_next(shashc_t* cursor)
{
    int bucket_index;
    SET_VALID_HASH;

    /* Cursor at end */
    if (cursor->bucket_entry == NULL)
        return 0;

    /* Still things in this bucket */
    if (bucket_next(cursor->bucket_entry)) {
        cursor->bucket_entry = bucket_next(cursor->bucket_entry);
        return 1;
    }

    /* Have to move to another bucket */
    bucket_index = cursor->bucket_number + 1;
    for (; bucket_index < hash->table_size; ++bucket_index) {
        if (hash->table[bucket_index] != NULL) {
            cursor->bucket_number = bucket_index;
            cursor->bucket_entry = hash->table[bucket_index];
            return 1;
        }
    }

    /* Otherwise, set to end/invalid entry */ 
    cursor->bucket_number = -1;
    cursor->bucket_entry = NULL;
    return 0;
}

/* Get the key of the current bucket.  NULL for end bucket.
 */
char* shashc_key(shashc_t* cursor)
{
    SET_VALID_HASH;
    if (cursor->bucket_entry)
        return bucket_key(cursor->bucket_entry);
    else
        return NULL;
}

/* Get the data of the current bucket.  NULL for end bucket.
 */
void* shashc_data(shashc_t* cursor)
{
    SET_VALID_HASH;
    if (cursor->bucket_entry)
        return bucket_value(cursor->bucket_entry);
    else
        return NULL;
}

/* True iff cursor is at the end.
 */
int shashc_end(shashc_t* cursor)
{
    SET_VALID_HASH;
    return (cursor->bucket_entry == NULL);
}
