/*source khash.h*/
/* The MIT License

   Copyright (c) 2008, 2009, 2011 by Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#ifndef __AC_KHASH_H
#define __AC_KHASH_H

/*!
  @header

  Generic hash table library.

  @copyright Heng Li
 */

#define AC_VERSION_KHASH_H "0.2.5"

#include <stdlib.h>
#include <string.h>
#include <limits.h>

/* compipler specific configuration */

#if UINT_MAX == 0xffffffffu
typedef unsigned int khint32_t;
#elif ULONG_MAX == 0xffffffffu
typedef unsigned long khint32_t;
#endif

#if ULONG_MAX == ULLONG_MAX
typedef unsigned long khint64_t;
#else
typedef unsigned long long khint64_t;
#endif

#ifdef _MSC_VER
#define inline __inline
#endif

typedef khint32_t khint_t;
typedef khint_t khiter_t;

#define __ac_HASH_PRIME_SIZE 32
static const khint32_t __ac_prime_list[__ac_HASH_PRIME_SIZE] =
{
  0ul,          3ul,          11ul,         23ul,         53ul,
  97ul,         193ul,        389ul,        769ul,        1543ul,
  3079ul,       6151ul,       12289ul,      24593ul,      49157ul,
  98317ul,      196613ul,     393241ul,     786433ul,     1572869ul,
  3145739ul,    6291469ul,    12582917ul,   25165843ul,   50331653ul,
  100663319ul,  201326611ul,  402653189ul,  805306457ul,  1610612741ul,
  3221225473ul, 4294967291ul
};

#define __ac_isempty(flag, i) ((flag[i>>4]>>((i&0xfU)<<1))&2)
#define __ac_isdel(flag, i) ((flag[i>>4]>>((i&0xfU)<<1))&1)
#define __ac_iseither(flag, i) ((flag[i>>4]>>((i&0xfU)<<1))&3)
#define __ac_set_isdel_false(flag, i) (flag[i>>4]&=~(1ul<<((i&0xfU)<<1)))
#define __ac_set_isempty_false(flag, i) (flag[i>>4]&=~(2ul<<((i&0xfU)<<1)))
#define __ac_set_isboth_false(flag, i) (flag[i>>4]&=~(3ul<<((i&0xfU)<<1)))
#define __ac_set_isdel_true(flag, i) (flag[i>>4]|=1ul<<((i&0xfU)<<1))

static const double __ac_HASH_UPPER = 0.77;

#define KHASH_DECLARE(name, khkey_t, khval_t)		 					\
	typedef struct {													\
		khint_t n_buckets, size, n_occupied, upper_bound;				\
		khint32_t *flags;												\
		khkey_t *keys;													\
		khval_t *vals;													\
	} kh_##name##_t;													\
	extern kh_##name##_t *kh_init_##name();								\
	extern void kh_destroy_##name(kh_##name##_t *h);					\
	extern void kh_clear_##name(kh_##name##_t *h);						\
	extern khint_t kh_get_##name(const kh_##name##_t *h, khkey_t key); 	\
	extern void kh_resize_##name(kh_##name##_t *h, khint_t new_n_buckets); \
	extern khint_t kh_put_##name(kh_##name##_t *h, khkey_t key, int *ret); \
	extern void kh_del_##name(kh_##name##_t *h, khint_t x);

#define KHASH_INIT2(name, SCOPE, khkey_t, khval_t, kh_is_map, __hash_func, __hash_equal) \
	typedef struct {													\
		khint_t n_buckets, size, n_occupied, upper_bound;				\
		khint32_t *flags;												\
		khkey_t *keys;													\
		khval_t *vals;													\
	} kh_##name##_t;													\
	SCOPE kh_##name##_t *kh_init_##name() {								\
		return (kh_##name##_t*)calloc(1, sizeof(kh_##name##_t));		\
	}																	\
	SCOPE void kh_destroy_##name(kh_##name##_t *h)						\
	{																	\
		if (h) {														\
			free(h->keys); free(h->flags);								\
			free(h->vals);												\
			free(h);													\
		}																\
	}																	\
	SCOPE void kh_clear_##name(kh_##name##_t *h)						\
	{																	\
		if (h && h->flags) {											\
			memset(h->flags, 0xaa, ((h->n_buckets>>4) + 1) * sizeof(khint32_t)); \
			h->size = h->n_occupied = 0;								\
		}																\
	}																	\
	SCOPE khint_t kh_get_##name(const kh_##name##_t *h, khkey_t key) 	\
	{																	\
		if (h->n_buckets) {												\
			khint_t inc, k, i, last;									\
			k = __hash_func(key); i = k % h->n_buckets;					\
			inc = 1 + k % (h->n_buckets - 1); last = i;					\
			while (!__ac_isempty(h->flags, i) && (__ac_isdel(h->flags, i) || !__hash_equal(h->keys[i], key))) { \
				if (i + inc >= h->n_buckets) i = i + inc - h->n_buckets; \
				else i += inc;											\
				if (i == last) return h->n_buckets;						\
			}															\
			return __ac_iseither(h->flags, i)? h->n_buckets : i;		\
		} else return 0;												\
	}																	\
	SCOPE void kh_resize_##name(kh_##name##_t *h, khint_t new_n_buckets) \
	{																	\
		khint32_t *new_flags = 0;										\
		khint_t j = 1;													\
		{																\
			khint_t t = __ac_HASH_PRIME_SIZE - 1;						\
			while (__ac_prime_list[t] > new_n_buckets) --t;				\
			new_n_buckets = __ac_prime_list[t+1];						\
			if (h->size >= (khint_t)(new_n_buckets * __ac_HASH_UPPER + 0.5)) j = 0;	\
			else {														\
				new_flags = (khint32_t*)malloc(((new_n_buckets>>4) + 1) * sizeof(khint32_t));	\
				memset(new_flags, 0xaa, ((new_n_buckets>>4) + 1) * sizeof(khint32_t)); \
				if (h->n_buckets < new_n_buckets) {						\
					h->keys = (khkey_t*)realloc(h->keys, new_n_buckets * sizeof(khkey_t)); \
					if (kh_is_map)										\
						h->vals = (khval_t*)realloc(h->vals, new_n_buckets * sizeof(khval_t)); \
				}														\
			}															\
		}																\
		if (j) {														\
			for (j = 0; j != h->n_buckets; ++j) {						\
				if (__ac_iseither(h->flags, j) == 0) {					\
					khkey_t key = h->keys[j];							\
					khval_t val;										\
					if (kh_is_map) val = h->vals[j];					\
					__ac_set_isdel_true(h->flags, j);					\
					while (1) {											\
						khint_t inc, k, i;								\
						k = __hash_func(key);							\
						i = k % new_n_buckets;							\
						inc = 1 + k % (new_n_buckets - 1);				\
						while (!__ac_isempty(new_flags, i)) {			\
							if (i + inc >= new_n_buckets) i = i + inc - new_n_buckets; \
							else i += inc;								\
						}												\
						__ac_set_isempty_false(new_flags, i);			\
						if (i < h->n_buckets && __ac_iseither(h->flags, i) == 0) { \
							{ khkey_t tmp = h->keys[i]; h->keys[i] = key; key = tmp; } \
							if (kh_is_map) { khval_t tmp = h->vals[i]; h->vals[i] = val; val = tmp; } \
							__ac_set_isdel_true(h->flags, i);			\
						} else {										\
							h->keys[i] = key;							\
							if (kh_is_map) h->vals[i] = val;			\
							break;										\
						}												\
					}													\
				}														\
			}															\
			if (h->n_buckets > new_n_buckets) {							\
				h->keys = (khkey_t*)realloc(h->keys, new_n_buckets * sizeof(khkey_t)); \
				if (kh_is_map)											\
					h->vals = (khval_t*)realloc(h->vals, new_n_buckets * sizeof(khval_t)); \
			}															\
			free(h->flags);												\
			h->flags = new_flags;										\
			h->n_buckets = new_n_buckets;								\
			h->n_occupied = h->size;									\
			h->upper_bound = (khint_t)(h->n_buckets * __ac_HASH_UPPER + 0.5); \
		}																\
	}																	\
	SCOPE khint_t kh_put_##name(kh_##name##_t *h, khkey_t key, int *ret) \
	{																	\
		khint_t x;														\
		if (h->n_occupied >= h->upper_bound) {							\
			if (h->n_buckets > (h->size<<1)) kh_resize_##name(h, h->n_buckets - 1); \
			else kh_resize_##name(h, h->n_buckets + 1);					\
		}																\
		{																\
			khint_t inc, k, i, site, last;								\
			x = site = h->n_buckets; k = __hash_func(key); i = k % h->n_buckets; \
			if (__ac_isempty(h->flags, i)) x = i;						\
			else {														\
				inc = 1 + k % (h->n_buckets - 1); last = i;				\
				while (!__ac_isempty(h->flags, i) && (__ac_isdel(h->flags, i) || !__hash_equal(h->keys[i], key))) { \
					if (__ac_isdel(h->flags, i)) site = i;				\
					if (i + inc >= h->n_buckets) i = i + inc - h->n_buckets; \
					else i += inc;										\
					if (i == last) { x = site; break; }					\
				}														\
				if (x == h->n_buckets) {								\
					if (__ac_isempty(h->flags, i) && site != h->n_buckets) x = site; \
					else x = i;											\
				}														\
			}															\
		}																\
		if (__ac_isempty(h->flags, x)) {								\
			h->keys[x] = key;											\
			__ac_set_isboth_false(h->flags, x);							\
			++h->size; ++h->n_occupied;									\
			*ret = 1;													\
		} else if (__ac_isdel(h->flags, x)) {							\
			h->keys[x] = key;											\
			__ac_set_isboth_false(h->flags, x);							\
			++h->size;													\
			*ret = 2;													\
		} else *ret = 0;												\
		return x;														\
	}																	\
	SCOPE void kh_del_##name(kh_##name##_t *h, khint_t x)				\
	{																	\
		if (x != h->n_buckets && !__ac_iseither(h->flags, x)) {			\
			__ac_set_isdel_true(h->flags, x);							\
			--h->size;													\
		}																\
	}

#define KHASH_INIT(name, khkey_t, khval_t, kh_is_map, __hash_func, __hash_equal) \
	KHASH_INIT2(name, static inline, khkey_t, khval_t, kh_is_map, __hash_func, __hash_equal)

/* --- BEGIN OF HASH FUNCTIONS --- */

/*! @function
  @abstract     Integer hash function
  @param  key   The integer [khint32_t]
  @return       The hash value [khint_t]
 */
#define kh_int_hash_func(key) (khint32_t)(key)
/*! @function
  @abstract     Integer comparison function
 */
#define kh_int_hash_equal(a, b) ((a) == (b))
/*! @function
  @abstract     64-bit integer hash function
  @param  key   The integer [khint64_t]
  @return       The hash value [khint_t]
 */
#define kh_int64_hash_func(key) (khint32_t)((key)>>33^(key)^(key)<<11)
/*! @function
  @abstract     64-bit integer comparison function
 */
#define kh_int64_hash_equal(a, b) ((a) == (b))
/*! @function
  @abstract     const char* hash function
  @param  s     Pointer to a null terminated string
  @return       The hash value
 */
static inline khint_t __ac_X31_hash_string(const char *s)
{
	khint_t h = *s;
	if (h) for (++s ; *s; ++s) h = (h << 5) - h + *s;
	return h;
}
/*! @function
  @abstract     Another interface to const char* hash function
  @param  key   Pointer to a null terminated string [const char*]
  @return       The hash value [khint_t]
 */
#define kh_str_hash_func(key) __ac_X31_hash_string(key)
/*! @function
  @abstract     Const char* comparison function
 */
#define kh_str_hash_equal(a, b) (strcmp(a, b) == 0)

/* --- END OF HASH FUNCTIONS --- */

/* Other necessary macros... */

/*!
  @abstract Type of the hash table.
  @param  name  Name of the hash table [symbol]
 */
#define khash_t(name) kh_##name##_t

/*! @function
  @abstract     Initiate a hash table.
  @param  name  Name of the hash table [symbol]
  @return       Pointer to the hash table [khash_t(name)*]
 */
#define kh_init(name) kh_init_##name()

/*! @function
  @abstract     Destroy a hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
 */
#define kh_destroy(name, h) kh_destroy_##name(h)

/*! @function
  @abstract     Reset a hash table without deallocating memory.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
 */
#define kh_clear(name, h) kh_clear_##name(h)

/*! @function
  @abstract     Resize a hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  s     New size [khint_t]
 */
#define kh_resize(name, h, s) kh_resize_##name(h, s)

/*! @function
  @abstract     Insert a key to the hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  k     Key [type of keys]
  @param  r     Extra return code: 0 if the key is present in the hash table;
                1 if the bucket is empty (never used); 2 if the element in
				the bucket has been deleted [int*]
  @return       Iterator to the inserted element [khint_t]
 */
#define kh_put(name, h, k, r) kh_put_##name(h, k, r)

/*! @function
  @abstract     Retrieve a key from the hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  k     Key [type of keys]
  @return       Iterator to the found element, or kh_end(h) is the element is absent [khint_t]
 */
#define kh_get(name, h, k) kh_get_##name(h, k)

/*! @function
  @abstract     Remove a key from the hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  k     Iterator to the element to be deleted [khint_t]
 */
#define kh_del(name, h, k) kh_del_##name(h, k)


/*! @function
  @abstract     Test whether a bucket contains data.
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  x     Iterator to the bucket [khint_t]
  @return       1 if containing data; 0 otherwise [int]
 */
#define kh_exist(h, x) (!__ac_iseither((h)->flags, (x)))

/*! @function
  @abstract     Get key given an iterator
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  x     Iterator to the bucket [khint_t]
  @return       Key [type of keys]
 */
#define kh_key(h, x) ((h)->keys[x])

/*! @function
  @abstract     Get value given an iterator
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  x     Iterator to the bucket [khint_t]
  @return       Value [type of values]
  @discussion   For hash sets, calling this results in segfault.
 */
#define kh_val(h, x) ((h)->vals[x])

/*! @function
  @abstract     Alias of kh_val()
 */
#define kh_value(h, x) ((h)->vals[x])

/*! @function
  @abstract     Get the start iterator
  @param  h     Pointer to the hash table [khash_t(name)*]
  @return       The start iterator [khint_t]
 */
#define kh_begin(h) (khint_t)(0)

/*! @function
  @abstract     Get the end iterator
  @param  h     Pointer to the hash table [khash_t(name)*]
  @return       The end iterator [khint_t]
 */
#define kh_end(h) ((h)->n_buckets)

/*! @function
  @abstract     Get the number of elements in the hash table
  @param  h     Pointer to the hash table [khash_t(name)*]
  @return       Number of elements in the hash table [khint_t]
 */
#define kh_size(h) ((h)->size)

/*! @function
  @abstract     Get the number of buckets in the hash table
  @param  h     Pointer to the hash table [khash_t(name)*]
  @return       Number of buckets in the hash table [khint_t]
 */
#define kh_n_buckets(h) ((h)->n_buckets)

/* More conenient interfaces */

/*! @function
  @abstract     Instantiate a hash set containing integer keys
  @param  name  Name of the hash table [symbol]
 */
#define KHASH_SET_INIT_INT(name)										\
	KHASH_INIT(name, khint32_t, char, 0, kh_int_hash_func, kh_int_hash_equal)

/*! @function
  @abstract     Instantiate a hash map containing integer keys
  @param  name  Name of the hash table [symbol]
  @param  khval_t  Type of values [type]
 */
#define KHASH_MAP_INIT_INT(name, khval_t)								\
	KHASH_INIT(name, khint32_t, khval_t, 1, kh_int_hash_func, kh_int_hash_equal)

/*! @function
  @abstract     Instantiate a hash map containing 64-bit integer keys
  @param  name  Name of the hash table [symbol]
 */
#define KHASH_SET_INIT_INT64(name)										\
	KHASH_INIT(name, khint64_t, char, 0, kh_int64_hash_func, kh_int64_hash_equal)

/*! @function
  @abstract     Instantiate a hash map containing 64-bit integer keys
  @param  name  Name of the hash table [symbol]
  @param  khval_t  Type of values [type]
 */
#define KHASH_MAP_INIT_INT64(name, khval_t)								\
	KHASH_INIT(name, khint64_t, khval_t, 1, kh_int64_hash_func, kh_int64_hash_equal)

typedef const char *kh_cstr_t;
/*! @function
  @abstract     Instantiate a hash map containing const char* keys
  @param  name  Name of the hash table [symbol]
 */
#define KHASH_SET_INIT_STR(name)										\
	KHASH_INIT(name, kh_cstr_t, char, 0, kh_str_hash_func, kh_str_hash_equal)

/*! @function
  @abstract     Instantiate a hash map containing const char* keys
  @param  name  Name of the hash table [symbol]
  @param  khval_t  Type of values [type]
 */
#define KHASH_MAP_INIT_STR(name, khval_t)								\
	KHASH_INIT(name, kh_cstr_t, khval_t, 1, kh_str_hash_func, kh_str_hash_equal)

#endif /* __AC_KHASH_H */


/*source knetfile.h and knetfile.c*/
/* The MIT License

   Copyright (c) 2008 by Genome Research Ltd (GRL).
                 2010 by Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Probably I will not do socket programming in the next few years and
   therefore I decide to heavily annotate this file, for Linux and
   Windows as well.  -ac */

#ifndef KNETFILE_H
#define KNETFILE_H

#include <stdint.h>
#include <fcntl.h>
#include <time.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <sys/types.h>

#ifndef _WIN32
#include <netdb.h>
#include <arpa/inet.h>
#include <sys/socket.h>
#define netread(fd, ptr, len) read(fd, ptr, len)
#define netwrite(fd, ptr, len) write(fd, ptr, len)
#define netclose(fd) close(fd)
#else
#include <winsock2.h>
#define netread(fd, ptr, len) recv(fd, ptr, len, 0)
#define netwrite(fd, ptr, len) send(fd, ptr, len, 0)
#define netclose(fd) closesocket(fd)
#endif

/* In winsock.h, the type of a socket is SOCKET, which is: "typedef
 * u_int SOCKET". An invalid SOCKET is: "(SOCKET)(~0)", or signed
 * integer -1. In knetfile.c, I use "int" for socket type
 * throughout. This should be improved to avoid confusion.
 *
 * In Linux/Mac, recv() and read() do almost the same thing. You can see
 * in the header file that netread() is simply an alias of read(). In
 * Windows, however, they are different and using recv() is mandatory.
 */

// FIXME: currently I/O is unbuffered

#define KNF_TYPE_LOCAL 1
#define KNF_TYPE_FTP   2
#define KNF_TYPE_HTTP  3

typedef struct knetFile_s {
  int type, fd;
  int64_t offset;
  char *host, *port;
  
  // the following are for FTP only
  int ctrl_fd, pasv_ip[4], pasv_port, max_response, no_reconnect, is_ready;
  char *response, *retr, *size_cmd;
	int64_t seek_offset; // for lazy seek
  int64_t file_size;

  // the following are for HTTP only
  char *path, *http_host;
} knetFile;

#define knet_tell(fp) ((fp)->offset)
#define knet_fileno(fp) ((fp)->fd)

#ifdef _WIN32
int knet_win32_init();
void knet_win32_destroy();
#endif

knetFile *knet_open(const char *fn, const char *mode);

/* 
   This only works with local files.
*/
knetFile *knet_dopen(int fd, const char *mode);

/*
  If ->is_ready==0, this routine updates ->fd; otherwise, it simply
  reads from ->fd.
*/
off_t knet_read(knetFile *fp, void *buf, off_t len);

/*
  This routine only sets ->offset and ->is_ready=0. It does not
  communicate with the FTP server.
*/
off_t knet_seek(knetFile *fp, int64_t off, int whence);
int knet_close(knetFile *fp);


/* This function tests if the file handler is ready for reading (or
 * writing if is_read==0). */
static int socket_wait(int fd, int is_read)
{
	fd_set fds, *fdr = 0, *fdw = 0;
	struct timeval tv;
	int ret;
	tv.tv_sec = 5; tv.tv_usec = 0; // 5 seconds time out
	FD_ZERO(&fds);
	FD_SET(fd, &fds);
	if (is_read) fdr = &fds;
	else fdw = &fds;
	ret = select(fd+1, fdr, fdw, 0, &tv);
#ifndef _WIN32
	if (ret == -1) perror("select");
#else
	if (ret == 0)
		fprintf(stderr, "select time-out\n");
	else if (ret == SOCKET_ERROR)
		fprintf(stderr, "select: %d\n", WSAGetLastError());
#endif
	return ret;
}

#ifndef _WIN32
/* This function does not work with Windows due to the lack of
 * getaddrinfo() in winsock. It is addapted from an example in "Beej's
 * Guide to Network Programming" (http://beej.us/guide/bgnet/). */
static int socket_connect(const char *host, const char *port)
{
#define __err_connect(func) do { perror(func); freeaddrinfo(res); return -1; } while (0)

	int on = 1, fd;
	struct linger lng = { 0, 0 };
	struct addrinfo hints, *res = 0;
	memset(&hints, 0, sizeof(struct addrinfo));
	hints.ai_family = AF_UNSPEC;
	hints.ai_socktype = SOCK_STREAM;
	/* In Unix/Mac, getaddrinfo() is the most convenient way to get
	 * server information. */
	if (getaddrinfo(host, port, &hints, &res) != 0) __err_connect("getaddrinfo");
	if ((fd = socket(res->ai_family, res->ai_socktype, res->ai_protocol)) == -1) __err_connect("socket");
	/* The following two setsockopt() are used by ftplib
	 * (http://nbpfaus.net/~pfau/ftplib/). I am not sure if they
	 * necessary. */
	if (setsockopt(fd, SOL_SOCKET, SO_REUSEADDR, &on, sizeof(on)) == -1) __err_connect("setsockopt");
	if (setsockopt(fd, SOL_SOCKET, SO_LINGER, &lng, sizeof(lng)) == -1) __err_connect("setsockopt");
	if (connect(fd, res->ai_addr, res->ai_addrlen) != 0) __err_connect("connect");
	freeaddrinfo(res);
	return fd;
}
#else
/* MinGW's printf has problem with "%lld" */
char *int64tostr(char *buf, int64_t x)
{
	int cnt;
	int i = 0;
	do {
		buf[i++] = '0' + x % 10;
		x /= 10;
	} while (x);
	buf[i] = 0;
	for (cnt = i, i = 0; i < cnt/2; ++i) {
		int c = buf[i]; buf[i] = buf[cnt-i-1]; buf[cnt-i-1] = c;
	}
	return buf;
}

int64_t strtoint64(const char *buf)
{
	int64_t x;
	for (x = 0; *buf != '\0'; ++buf)
		x = x * 10 + ((int64_t) *buf - 48);
	return x;
}
/* In windows, the first thing is to establish the TCP connection. */
int knet_win32_init()
{
	WSADATA wsaData;
	return WSAStartup(MAKEWORD(2, 2), &wsaData);
}
void knet_win32_destroy()
{
	WSACleanup();
}
/* A slightly modfied version of the following function also works on
 * Mac (and presummably Linux). However, this function is not stable on
 * my Mac. It sometimes works fine but sometimes does not. Therefore for
 * non-Windows OS, I do not use this one. */
static SOCKET socket_connect(const char *host, const char *port)
{
#define __err_connect(func)										\
	do {														\
		fprintf(stderr, "%s: %d\n", func, WSAGetLastError());	\
		return -1;												\
	} while (0)

	int on = 1;
	SOCKET fd;
	struct linger lng = { 0, 0 };
	struct sockaddr_in server;
	struct hostent *hp = 0;
	// open socket
	if ((fd = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP)) == INVALID_SOCKET) __err_connect("socket");
	if (setsockopt(fd, SOL_SOCKET, SO_REUSEADDR, (char*)&on, sizeof(on)) == -1) __err_connect("setsockopt");
	if (setsockopt(fd, SOL_SOCKET, SO_LINGER, (char*)&lng, sizeof(lng)) == -1) __err_connect("setsockopt");
	// get host info
	if (isalpha(host[0])) hp = gethostbyname(host);
	else {
		struct in_addr addr;
		addr.s_addr = inet_addr(host);
		hp = gethostbyaddr((char*)&addr, 4, AF_INET);
	}
	if (hp == 0) __err_connect("gethost");
	// connect
	server.sin_addr.s_addr = *((unsigned long*)hp->h_addr);
	server.sin_family= AF_INET;
	server.sin_port = htons(atoi(port));
	if (connect(fd, (struct sockaddr*)&server, sizeof(server)) != 0) __err_connect("connect");
	// freehostent(hp); // strangely in MSDN, hp is NOT freed (memory leak?!)
	return fd;
}
#endif

static off_t my_netread(int fd, void *buf, off_t len)
{
	off_t rest = len, curr, l = 0;
	/* recv() and read() may not read the required length of data with
	 * one call. They have to be called repeatedly. */
	while (rest) {
		if (socket_wait(fd, 1) <= 0) break; // socket is not ready for reading
		curr = netread(fd, buf + l, rest);
		/* According to the glibc manual, section 13.2, a zero returned
		 * value indicates end-of-file (EOF), which should mean that
		 * read() will not return zero if EOF has not been met but data
		 * are not immediately available. */
		if (curr == 0) break;
		l += curr; rest -= curr;
	}
	return l;
}

/*************************
 * FTP specific routines *
 *************************/

static int kftp_get_response(knetFile *ftp)
{
#ifndef _WIN32
	unsigned char c;
#else
	char c;
#endif
	int n = 0;
	char *p;
	if (socket_wait(ftp->ctrl_fd, 1) <= 0) return 0;
	while (netread(ftp->ctrl_fd, &c, 1)) { // FIXME: this is *VERY BAD* for unbuffered I/O
		//fputc(c, stderr);
		if (n >= ftp->max_response) {
			ftp->max_response = ftp->max_response? ftp->max_response<<1 : 256;
			ftp->response = realloc(ftp->response, ftp->max_response);
		}
		ftp->response[n++] = c;
		if (c == '\n') {
			if (n >= 4 && isdigit(ftp->response[0]) && isdigit(ftp->response[1]) && isdigit(ftp->response[2])
				&& ftp->response[3] != '-') break;
			n = 0;
			continue;
		}
	}
	if (n < 2) return -1;
	ftp->response[n-2] = 0;
	return strtol(ftp->response, &p, 0);
}

static int kftp_send_cmd(knetFile *ftp, const char *cmd, int is_get)
{
	if (socket_wait(ftp->ctrl_fd, 0) <= 0) return -1; // socket is not ready for writing
	netwrite(ftp->ctrl_fd, cmd, strlen(cmd));
	return is_get? kftp_get_response(ftp) : 0;
}

static int kftp_pasv_prep(knetFile *ftp)
{
	char *p;
	int v[6];
	kftp_send_cmd(ftp, "PASV\r\n", 1);
	for (p = ftp->response; *p && *p != '('; ++p);
	if (*p != '(') return -1;
	++p;
	sscanf(p, "%d,%d,%d,%d,%d,%d", &v[0], &v[1], &v[2], &v[3], &v[4], &v[5]);
	memcpy(ftp->pasv_ip, v, 4 * sizeof(int));
	ftp->pasv_port = (v[4]<<8&0xff00) + v[5];
	return 0;
}


static int kftp_pasv_connect(knetFile *ftp)
{
	char host[80], port[10];
	if (ftp->pasv_port == 0) {
		fprintf(stderr, "[kftp_pasv_connect] kftp_pasv_prep() is not called before hand.\n");
		return -1;
	}
	sprintf(host, "%d.%d.%d.%d", ftp->pasv_ip[0], ftp->pasv_ip[1], ftp->pasv_ip[2], ftp->pasv_ip[3]);
	sprintf(port, "%d", ftp->pasv_port);
	ftp->fd = socket_connect(host, port);
	if (ftp->fd == -1) return -1;
	return 0;
}

int kftp_connect(knetFile *ftp)
{
	ftp->ctrl_fd = socket_connect(ftp->host, ftp->port);
	if (ftp->ctrl_fd == -1) return -1;
	kftp_get_response(ftp);
	kftp_send_cmd(ftp, "USER anonymous\r\n", 1);
	kftp_send_cmd(ftp, "PASS kftp@\r\n", 1);
	kftp_send_cmd(ftp, "TYPE I\r\n", 1);
	return 0;
}

int kftp_reconnect(knetFile *ftp)
{
	if (ftp->ctrl_fd != -1) {
		netclose(ftp->ctrl_fd);
		ftp->ctrl_fd = -1;
	}
	netclose(ftp->fd);
	ftp->fd = -1;
	return kftp_connect(ftp);
}

// initialize ->type, ->host, ->retr and ->size
knetFile *kftp_parse_url(const char *fn, const char *mode)
{
	knetFile *fp;
	char *p;
	int l;
	if (strstr(fn, "ftp://") != fn) return 0;
	for (p = (char*)fn + 6; *p && *p != '/'; ++p);
	if (*p != '/') return 0;
	l = p - fn - 6;
	fp = calloc(1, sizeof(knetFile));
	fp->type = KNF_TYPE_FTP;
	fp->fd = -1;
	/* the Linux/Mac version of socket_connect() also recognizes a port
	 * like "ftp", but the Windows version does not. */
	fp->port = strdup("21");
	fp->host = calloc(l + 1, 1);
	if (strchr(mode, 'c')) fp->no_reconnect = 1;
	strncpy(fp->host, fn + 6, l);
	fp->retr = calloc(strlen(p) + 8, 1);
	sprintf(fp->retr, "RETR %s\r\n", p);
    fp->size_cmd = calloc(strlen(p) + 8, 1);
    sprintf(fp->size_cmd, "SIZE %s\r\n", p);
	fp->seek_offset = 0;
	return fp;
}
// place ->fd at offset off
int kftp_connect_file(knetFile *fp)
{
	int ret;
	long long file_size;
	if (fp->fd != -1) {
		netclose(fp->fd);
		if (fp->no_reconnect) kftp_get_response(fp);
	}
	kftp_pasv_prep(fp);
    kftp_send_cmd(fp, fp->size_cmd, 1);
#ifndef _WIN32
    if ( sscanf(fp->response,"%*d %lld", &file_size) != 1 )
    {
        fprintf(stderr,"[kftp_connect_file] %s\n", fp->response);
        return -1;
    }
#else
	const char *p = fp->response;
	while (*p != ' ') ++p;
	while (*p < '0' || *p > '9') ++p;
	file_size = strtoint64(p);
#endif
	fp->file_size = file_size;
	if (fp->offset>=0) {
		char tmp[32];
#ifndef _WIN32
		sprintf(tmp, "REST %lld\r\n", (long long)fp->offset);
#else
		strcpy(tmp, "REST ");
		int64tostr(tmp + 5, fp->offset);
		strcat(tmp, "\r\n");
#endif
		kftp_send_cmd(fp, tmp, 1);
	}
	kftp_send_cmd(fp, fp->retr, 0);
	kftp_pasv_connect(fp);
	ret = kftp_get_response(fp);
	if (ret != 150) {
		fprintf(stderr, "[kftp_connect_file] %s\n", fp->response);
		netclose(fp->fd);
		fp->fd = -1;
		return -1;
	}
	fp->is_ready = 1;
	return 0;
}


/**************************
 * HTTP specific routines *
 **************************/

knetFile *khttp_parse_url(const char *fn, const char *mode)
{
	knetFile *fp;
	char *p, *proxy, *q;
	int l;
	if (strstr(fn, "http://") != fn) return 0;
	// set ->http_host
	for (p = (char*)fn + 7; *p && *p != '/'; ++p);
	l = p - fn - 7;
	fp = calloc(1, sizeof(knetFile));
	fp->http_host = calloc(l + 1, 1);
	strncpy(fp->http_host, fn + 7, l);
	fp->http_host[l] = 0;
	for (q = fp->http_host; *q && *q != ':'; ++q);
	if (*q == ':') *q++ = 0;
	// get http_proxy
	proxy = getenv("http_proxy");
	// set ->host, ->port and ->path
	if (proxy == 0) {
		fp->host = strdup(fp->http_host); // when there is no proxy, server name is identical to http_host name.
		fp->port = strdup(*q? q : "80");
		fp->path = strdup(*p? p : "/");
	} else {
		fp->host = (strstr(proxy, "http://") == proxy)? strdup(proxy + 7) : strdup(proxy);
		for (q = fp->host; *q && *q != ':'; ++q);
		if (*q == ':') *q++ = 0; 
		fp->port = strdup(*q? q : "80");
		fp->path = strdup(fn);
	}
	fp->type = KNF_TYPE_HTTP;
	fp->ctrl_fd = fp->fd = -1;
	fp->seek_offset = 0;
	return fp;
}

int khttp_connect_file(knetFile *fp)
{
	int ret, l = 0;
	char *buf, *p;
	if (fp->fd != -1) netclose(fp->fd);
	fp->fd = socket_connect(fp->host, fp->port);
	buf = calloc(0x10000, 1); // FIXME: I am lazy... But in principle, 64KB should be large enough.
	l += sprintf(buf + l, "GET %s HTTP/1.0\r\nHost: %s\r\n", fp->path, fp->http_host);
    l += sprintf(buf + l, "Range: bytes=%lld-\r\n", (long long)fp->offset);
	l += sprintf(buf + l, "\r\n");
	netwrite(fp->fd, buf, l);
	l = 0;
	while (netread(fp->fd, buf + l, 1)) { // read HTTP header; FIXME: bad efficiency
		if (buf[l] == '\n' && l >= 3)
			if (strncmp(buf + l - 3, "\r\n\r\n", 4) == 0) break;
		++l;
	}
	buf[l] = 0;
	if (l < 14) { // prematured header
		netclose(fp->fd);
		fp->fd = -1;
		return -1;
	}
	ret = strtol(buf + 8, &p, 0); // HTTP return code
	if (ret == 200 && fp->offset>0) { // 200 (complete result); then skip beginning of the file
		off_t rest = fp->offset;
		while (rest) {
			off_t l = rest < 0x10000? rest : 0x10000;
			rest -= my_netread(fp->fd, buf, l);
		}
	} else if (ret != 206 && ret != 200) {
		free(buf);
		fprintf(stderr, "[khttp_connect_file] fail to open file (HTTP code: %d).\n", ret);
		netclose(fp->fd);
		fp->fd = -1;
		return -1;
	}
	free(buf);
	fp->is_ready = 1;
	return 0;
}

/********************
 * Generic routines *
 ********************/

knetFile *knet_open(const char *fn, const char *mode)
{
	knetFile *fp = 0;
	if (mode[0] != 'r') {
		fprintf(stderr, "[kftp_open] only mode \"r\" is supported.\n");
		return 0;
	}
	if (strstr(fn, "ftp://") == fn) {
		fp = kftp_parse_url(fn, mode);
		if (fp == 0) return 0;
		if (kftp_connect(fp) == -1) {
			knet_close(fp);
			return 0;
		}
		kftp_connect_file(fp);
	} else if (strstr(fn, "http://") == fn) {
		fp = khttp_parse_url(fn, mode);
		if (fp == 0) return 0;
		khttp_connect_file(fp);
	} else { // local file
#ifdef _WIN32
		/* In windows, O_BINARY is necessary. In Linux/Mac, O_BINARY may
		 * be undefined on some systems, although it is defined on my
		 * Mac and the Linux I have tested on. */
		int fd = open(fn, O_RDONLY | O_BINARY);
#else		
		int fd = open(fn, O_RDONLY);
#endif
		if (fd == -1) {
			perror("open");
			return 0;
		}
		fp = (knetFile*)calloc(1, sizeof(knetFile));
		fp->type = KNF_TYPE_LOCAL;
		fp->fd = fd;
		fp->ctrl_fd = -1;
	}
	if (fp && fp->fd == -1) {
		knet_close(fp);
		return 0;
	}
	return fp;
}

knetFile *knet_dopen(int fd, const char *mode)
{
	knetFile *fp = (knetFile*)calloc(1, sizeof(knetFile));
	fp->type = KNF_TYPE_LOCAL;
	fp->fd = fd;
	return fp;
}

off_t knet_read(knetFile *fp, void *buf, off_t len)
{
	off_t l = 0;
	if (fp->fd == -1) return 0;
	if (fp->type == KNF_TYPE_FTP) {
		if (fp->is_ready == 0) {
			if (!fp->no_reconnect) kftp_reconnect(fp);
			kftp_connect_file(fp);
		}
	} else if (fp->type == KNF_TYPE_HTTP) {
		if (fp->is_ready == 0)
			khttp_connect_file(fp);
	}
	if (fp->type == KNF_TYPE_LOCAL) { // on Windows, the following block is necessary; not on UNIX
		off_t rest = len, curr;
		while (rest) {
			do {
				curr = read(fp->fd, buf + l, rest);
			} while (curr < 0 && EINTR == errno);
			if (curr < 0) return -1;
			if (curr == 0) break;
			l += curr; rest -= curr;
		}
	} else l = my_netread(fp->fd, buf, len);
	fp->offset += l;
	return l;
}

off_t knet_seek(knetFile *fp, int64_t off, int whence)
{
	if (whence == SEEK_SET && off == fp->offset) return 0;
	if (fp->type == KNF_TYPE_LOCAL) {
		/* Be aware that lseek() returns the offset after seeking,
		 * while fseek() returns zero on success. */
		off_t offset = lseek(fp->fd, off, whence);
		if (offset == -1) {
            // Be silent, it is OK for knet_seek to fail when the file is streamed
            // fprintf(stderr,"[knet_seek] %s\n", strerror(errno));
			return -1;
		}
		fp->offset = offset;
		return 0;
	}
    else if (fp->type == KNF_TYPE_FTP) 
    {
        if (whence==SEEK_CUR)
            fp->offset += off;
        else if (whence==SEEK_SET)
            fp->offset = off;
        else if ( whence==SEEK_END)
            fp->offset = fp->file_size+off;
		fp->is_ready = 0;
		return 0;
	} 
    else if (fp->type == KNF_TYPE_HTTP) 
    {
		if (whence == SEEK_END) { // FIXME: can we allow SEEK_END in future?
			fprintf(stderr, "[knet_seek] SEEK_END is not supported for HTTP. Offset is unchanged.\n");
			errno = ESPIPE;
			return -1;
		}
        if (whence==SEEK_CUR)
            fp->offset += off;
        else if (whence==SEEK_SET)
            fp->offset = off;
		fp->is_ready = 0;
		return 0;
	}
	errno = EINVAL;
    fprintf(stderr,"[knet_seek] %s\n", strerror(errno));
	return -1;
}

int knet_close(knetFile *fp)
{
	if (fp == 0) return 0;
	if (fp->ctrl_fd != -1) netclose(fp->ctrl_fd); // FTP specific
	if (fp->fd != -1) {
		/* On Linux/Mac, netclose() is an alias of close(), but on
		 * Windows, it is an alias of closesocket(). */
		if (fp->type == KNF_TYPE_LOCAL) close(fp->fd);
		else netclose(fp->fd);
	}
	free(fp->host); free(fp->port);
	free(fp->response); free(fp->retr); // FTP specific
	free(fp->path); free(fp->http_host); // HTTP specific
	free(fp);
	return 0;
}

#endif /* KNETFILE_H */


/*source razf.h and rasf.c*/
 /*-
 * RAZF : Random Access compressed(Z) File
 * Version: 1.0
 * Release Date: 2008-10-27
 *
 * Copyright 2008, Jue Ruan <ruanjue@gmail.com>, Heng Li <lh3@sanger.ac.uk>
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */


#ifndef __RAZF_RJ_H
#define __RAZF_RJ_H

#include <stdint.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include "zlib.h"

#if ZLIB_VERNUM < 0x1221
#define _RZ_READONLY
struct _gz_header_s;
typedef struct _gz_header_s _gz_header;
#define gz_header _gz_header
#endif

#define WINDOW_BITS   15

#ifndef RZ_BLOCK_SIZE
#define RZ_BLOCK_SIZE (1<<WINDOW_BITS)
#endif

#ifndef RZ_BUFFER_SIZE
#define RZ_BUFFER_SIZE 4096
#endif

#ifndef RZ_COMPRESS_LEVEL
#define RZ_COMPRESS_LEVEL 6
#endif

#define RZ_BIN_SIZE ((1LLU << 32) / RZ_BLOCK_SIZE)

typedef struct {
	uint32_t *cell_offsets; // i
	int64_t  *bin_offsets; // i / BIN_SIZE
	int size;
	int cap;
} ZBlockIndex;
/* When storing index, output bytes in Big-Endian everywhere */

#define FILE_TYPE_RZ	1
#define FILE_TYPE_PLAIN	2
#define FILE_TYPE_GZ	3

typedef struct RandomAccessZFile  {
	char mode; /* 'w' : write mode; 'r' : read mode */
	int file_type;
	/* plain file or rz file, razf_read support plain file as input too, in this case, razf_read work as buffered fread */
#ifdef _USE_KNETFILE
    union {
        knetFile *fpr;
        int fpw;
    } x;
#else
	int filedes; /* the file descriptor */
#endif
	z_stream *stream;
	ZBlockIndex *index;
	int64_t in, out, end, src_end;
	/* in: n bytes total in; out: n bytes total out; */
	/* end: the end of all data blocks, while the start of index; src_end: the true end position in uncompressed file */
	int buf_flush; // buffer should be flush, suspend inflate util buffer is empty
	int64_t block_pos, block_off, next_block_pos;
	/* block_pos: the start postiion of current block  in compressed file */
	/* block_off: tell how many bytes have been read from current block */
	void *inbuf, *outbuf;
	int header_size;
	gz_header *header;
	/* header is used to transfer inflate_state->mode from HEAD to TYPE after call inflateReset */
	int buf_off, buf_len;
	int z_err, z_eof;
	int seekable;
	/* Indice where the source is seekable */
	int load_index;
	/* set has_index to 0 in mode 'w', then index will be discarded */
} RAZF;


RAZF* razf_dopen(int data_fd, const char *mode);
RAZF *razf_open(const char *fn, const char *mode);
int razf_write(RAZF* rz, const void *data, int size);
int razf_read(RAZF* rz, void *data, int size);
int64_t razf_seek(RAZF* rz, int64_t pos, int where);
void razf_close(RAZF* rz);

#define razf_tell(rz) ((rz)->out)

RAZF* razf_open2(const char *filename, const char *mode);
RAZF* razf_dopen2(int fd, const char *mode);
uint64_t razf_tell2(RAZF *rz);
int64_t razf_seek2(RAZF *rz, uint64_t voffset, int where);

#if ZLIB_VERNUM < 0x1221
struct _gz_header_s {
    int     text;
    uLong   time;
    int     xflags;
    int     os;
    Bytef   *extra;
    uInt    extra_len;
    uInt    extra_max;
    Bytef   *name;
    uInt    name_max;
    Bytef   *comment;
    uInt    comm_max;
    int     hcrc;
    int     done;
};
#warning "zlib < 1.2.2.1; RAZF writing is disabled."
#endif

#define DEF_MEM_LEVEL 8

static inline uint32_t byte_swap_4(uint32_t v){
	v = ((v & 0x0000FFFFU) << 16) | (v >> 16);
	return ((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8);
}

static inline uint64_t byte_swap_8(uint64_t v){
	v = ((v & 0x00000000FFFFFFFFLLU) << 32) | (v >> 32);
	v = ((v & 0x0000FFFF0000FFFFLLU) << 16) | ((v & 0xFFFF0000FFFF0000LLU) >> 16);
	return ((v & 0x00FF00FF00FF00FFLLU) << 8) | ((v & 0xFF00FF00FF00FF00LLU) >> 8);
}

static inline int is_big_endian(){
	int x = 0x01;
	char *c = (char*)&x;
	return (c[0] != 0x01);
}

#ifndef _RZ_READONLY
static void add_zindex(RAZF *rz, int64_t in, int64_t out){
	if(rz->index->size == rz->index->cap){
		rz->index->cap = rz->index->cap * 1.5 + 2;
		rz->index->cell_offsets = realloc(rz->index->cell_offsets, sizeof(int) * rz->index->cap);
		rz->index->bin_offsets  = realloc(rz->index->bin_offsets, sizeof(int64_t) * (rz->index->cap/RZ_BIN_SIZE + 1));
	}
	if(rz->index->size % RZ_BIN_SIZE == 0) rz->index->bin_offsets[rz->index->size / RZ_BIN_SIZE] = out;
	rz->index->cell_offsets[rz->index->size] = out - rz->index->bin_offsets[rz->index->size / RZ_BIN_SIZE];
	rz->index->size ++;
}

static void save_zindex(RAZF *rz, int fd){
	int32_t i, v32;
	int is_be;
	is_be = is_big_endian();
	if(is_be) write(fd, &rz->index->size, sizeof(int));
	else {
		v32 = byte_swap_4((uint32_t)rz->index->size);
		write(fd, &v32, sizeof(uint32_t));
	}
	v32 = rz->index->size / RZ_BIN_SIZE + 1;
	if(!is_be){
		for(i=0;i<v32;i++) rz->index->bin_offsets[i]  = byte_swap_8((uint64_t)rz->index->bin_offsets[i]);
		for(i=0;i<rz->index->size;i++) rz->index->cell_offsets[i] = byte_swap_4((uint32_t)rz->index->cell_offsets[i]);
	}
	write(fd, rz->index->bin_offsets, sizeof(int64_t) * v32);
	write(fd, rz->index->cell_offsets, sizeof(int32_t) * rz->index->size);
}
#endif

#ifdef _USE_KNETFILE
static void load_zindex(RAZF *rz, knetFile *fp){
#else
static void load_zindex(RAZF *rz, int fd){
#endif
	int32_t i, v32;
	int is_be;
	if(!rz->load_index) return;
	if(rz->index == NULL) rz->index = malloc(sizeof(ZBlockIndex));
	is_be = is_big_endian();
#ifdef _USE_KNETFILE
	knet_read(fp, &rz->index->size, sizeof(int));
#else
	read(fd, &rz->index->size, sizeof(int));
#endif
	if(!is_be) rz->index->size = byte_swap_4((uint32_t)rz->index->size);
	rz->index->cap = rz->index->size;
	v32 = rz->index->size / RZ_BIN_SIZE + 1;
	rz->index->bin_offsets  = malloc(sizeof(int64_t) * v32);
#ifdef _USE_KNETFILE
	knet_read(fp, rz->index->bin_offsets, sizeof(int64_t) * v32);
#else
	read(fd, rz->index->bin_offsets, sizeof(int64_t) * v32);
#endif
	rz->index->cell_offsets = malloc(sizeof(int) * rz->index->size);
#ifdef _USE_KNETFILE
	knet_read(fp, rz->index->cell_offsets, sizeof(int) * rz->index->size);
#else
	read(fd, rz->index->cell_offsets, sizeof(int) * rz->index->size);
#endif
	if(!is_be){
		for(i=0;i<v32;i++) rz->index->bin_offsets[i] = byte_swap_8((uint64_t)rz->index->bin_offsets[i]);
		for(i=0;i<rz->index->size;i++) rz->index->cell_offsets[i] = byte_swap_4((uint32_t)rz->index->cell_offsets[i]);
	}
}

#ifdef _RZ_READONLY
static RAZF* razf_open_w(int fd)
{
	fprintf(stderr, "[razf_open_w] Writing is not available with zlib ver < 1.2.2.1\n");
	return 0;
}
#else
static RAZF* razf_open_w(int fd){
	RAZF *rz;
#ifdef _WIN32
	setmode(fd, O_BINARY);
#endif
	rz = calloc(1, sizeof(RAZF));
	rz->mode = 'w';
#ifdef _USE_KNETFILE
    rz->x.fpw = fd;
#else
	rz->filedes = fd;
#endif
	rz->stream = calloc(sizeof(z_stream), 1);
	rz->inbuf  = malloc(RZ_BUFFER_SIZE);
	rz->outbuf = malloc(RZ_BUFFER_SIZE);
	rz->index = calloc(sizeof(ZBlockIndex), 1);
	deflateInit2(rz->stream, RZ_COMPRESS_LEVEL, Z_DEFLATED, WINDOW_BITS + 16, DEF_MEM_LEVEL, Z_DEFAULT_STRATEGY);
	rz->stream->avail_out = RZ_BUFFER_SIZE;
	rz->stream->next_out  = rz->outbuf;
	rz->header = calloc(sizeof(gz_header), 1);
	rz->header->os    = 0x03; //Unix
	rz->header->text  = 0;
	rz->header->time  = 0;
	rz->header->extra = malloc(7);
	strncpy((char*)rz->header->extra, "RAZF", 4);
	rz->header->extra[4] = 1; // obsolete field
	// block size = RZ_BLOCK_SIZE, Big-Endian
	rz->header->extra[5] = RZ_BLOCK_SIZE >> 8;
	rz->header->extra[6] = RZ_BLOCK_SIZE & 0xFF;
	rz->header->extra_len = 7;
	rz->header->name = rz->header->comment  = 0;
	rz->header->hcrc = 0;
	deflateSetHeader(rz->stream, rz->header);
	rz->block_pos = rz->block_off = 0;
	return rz;
}

static void _razf_write(RAZF* rz, const void *data, int size){
	int tout;
	rz->stream->avail_in = size;
	rz->stream->next_in  = (void*)data;
	while(1){
		tout = rz->stream->avail_out;
		deflate(rz->stream, Z_NO_FLUSH);
		rz->out += tout - rz->stream->avail_out;
		if(rz->stream->avail_out) break;
#ifdef _USE_KNETFILE
		write(rz->x.fpw, rz->outbuf, RZ_BUFFER_SIZE - rz->stream->avail_out);
#else
		write(rz->filedes, rz->outbuf, RZ_BUFFER_SIZE - rz->stream->avail_out);
#endif
		rz->stream->avail_out = RZ_BUFFER_SIZE;
		rz->stream->next_out  = rz->outbuf;
		if(rz->stream->avail_in == 0) break;
	};
	rz->in += size - rz->stream->avail_in;
	rz->block_off += size - rz->stream->avail_in;
}

static void razf_flush(RAZF *rz){
	uint32_t tout;
	if(rz->buf_len){
		_razf_write(rz, rz->inbuf, rz->buf_len);
		rz->buf_off = rz->buf_len = 0;
	}
	if(rz->stream->avail_out){
#ifdef _USE_KNETFILE    
		write(rz->x.fpw, rz->outbuf, RZ_BUFFER_SIZE - rz->stream->avail_out);
#else        
		write(rz->filedes, rz->outbuf, RZ_BUFFER_SIZE - rz->stream->avail_out);
#endif
		rz->stream->avail_out = RZ_BUFFER_SIZE;
		rz->stream->next_out  = rz->outbuf;
	}
	while(1){
		tout = rz->stream->avail_out;
		deflate(rz->stream, Z_FULL_FLUSH);
		rz->out += tout - rz->stream->avail_out;
		if(rz->stream->avail_out == 0){
#ifdef _USE_KNETFILE    
			write(rz->x.fpw, rz->outbuf, RZ_BUFFER_SIZE - rz->stream->avail_out);
#else            
			write(rz->filedes, rz->outbuf, RZ_BUFFER_SIZE - rz->stream->avail_out);
#endif
			rz->stream->avail_out = RZ_BUFFER_SIZE;
			rz->stream->next_out  = rz->outbuf;
		} else break;
	}
	rz->block_pos = rz->out;
	rz->block_off = 0;
}

static void razf_end_flush(RAZF *rz){
	uint32_t tout;
	if(rz->buf_len){
		_razf_write(rz, rz->inbuf, rz->buf_len);
		rz->buf_off = rz->buf_len = 0;
	}
	while(1){
		tout = rz->stream->avail_out;
		deflate(rz->stream, Z_FINISH);
		rz->out += tout - rz->stream->avail_out;
		if(rz->stream->avail_out < RZ_BUFFER_SIZE){
#ifdef _USE_KNETFILE        
			write(rz->x.fpw, rz->outbuf, RZ_BUFFER_SIZE - rz->stream->avail_out);
#else            
			write(rz->filedes, rz->outbuf, RZ_BUFFER_SIZE - rz->stream->avail_out);
#endif
			rz->stream->avail_out = RZ_BUFFER_SIZE;
			rz->stream->next_out  = rz->outbuf;
		} else break;
	}
}

static void _razf_buffered_write(RAZF *rz, const void *data, int size){
	int i, n;
	while(1){
		if(rz->buf_len == RZ_BUFFER_SIZE){
			_razf_write(rz, rz->inbuf, rz->buf_len);
			rz->buf_len = 0;
		}
		if(size + rz->buf_len < RZ_BUFFER_SIZE){
			for(i=0;i<size;i++) ((char*)rz->inbuf + rz->buf_len)[i] = ((char*)data)[i];
			rz->buf_len += size;
			return;
		} else {
			n = RZ_BUFFER_SIZE - rz->buf_len;
			for(i=0;i<n;i++) ((char*)rz->inbuf + rz->buf_len)[i] = ((char*)data)[i];
			size -= n;
			data += n;
			rz->buf_len += n;
		}
	}
}

int razf_write(RAZF* rz, const void *data, int size){
	int ori_size, n;
	int64_t next_block;
	ori_size = size;
	next_block = ((rz->in / RZ_BLOCK_SIZE) + 1) * RZ_BLOCK_SIZE;
	while(rz->in + rz->buf_len + size >= next_block){
		n = next_block - rz->in - rz->buf_len;
		_razf_buffered_write(rz, data, n);
		data += n;
		size -= n;
		razf_flush(rz);
		add_zindex(rz, rz->in, rz->out);
		next_block = ((rz->in / RZ_BLOCK_SIZE) + 1) * RZ_BLOCK_SIZE;
	}
	_razf_buffered_write(rz, data, size);
	return ori_size;
}
#endif

/* gzip flag byte */
#define ASCII_FLAG   0x01 /* bit 0 set: file probably ascii text */
#define HEAD_CRC     0x02 /* bit 1 set: header CRC present */
#define EXTRA_FIELD  0x04 /* bit 2 set: extra field present */
#define ORIG_NAME    0x08 /* bit 3 set: original file name present */
#define COMMENT      0x10 /* bit 4 set: file comment present */
#define RESERVED     0xE0 /* bits 5..7: reserved */

static int _read_gz_header(unsigned char *data, int size, int *extra_off, int *extra_len){
	int method, flags, n, len;
	if(size < 2) return 0;
	if(data[0] != 0x1f || data[1] != 0x8b) return 0;
	if(size < 4) return 0;
	method = data[2];
	flags  = data[3];
	if(method != Z_DEFLATED || (flags & RESERVED)) return 0;
	n = 4 + 6; // Skip 6 bytes
	*extra_off = n + 2;
	*extra_len = 0;
	if(flags & EXTRA_FIELD){
		if(size < n + 2) return 0;
		len = ((int)data[n + 1] << 8) | data[n];
		n += 2;
		*extra_off = n;
		while(len){
			if(n >= size) return 0;
			n ++;
			len --;
		}
		*extra_len = n - (*extra_off);
	}
	if(flags & ORIG_NAME) while(n < size && data[n++]);
	if(flags & COMMENT) while(n < size && data[n++]);
	if(flags & HEAD_CRC){
		if(n + 2 > size) return 0;
		n += 2;
	}
	return n;
}

#ifdef _USE_KNETFILE
static RAZF* razf_open_r(knetFile *fp, int _load_index){
#else
static RAZF* razf_open_r(int fd, int _load_index){
#endif
	RAZF *rz;
	int ext_off, ext_len;
	int n, is_be, ret;
	int64_t end;
	unsigned char c[] = "RAZF";
	rz = calloc(1, sizeof(RAZF));
	rz->mode = 'r';
#ifdef _USE_KNETFILE
    rz->x.fpr = fp;
#else
#ifdef _WIN32
	setmode(fd, O_BINARY);
#endif
	rz->filedes = fd;
#endif
	rz->stream = calloc(sizeof(z_stream), 1);
	rz->inbuf  = malloc(RZ_BUFFER_SIZE);
	rz->outbuf = malloc(RZ_BUFFER_SIZE);
	rz->end = rz->src_end = 0x7FFFFFFFFFFFFFFFLL;
#ifdef _USE_KNETFILE
    n = knet_read(rz->x.fpr, rz->inbuf, RZ_BUFFER_SIZE);
#else
	n = read(rz->filedes, rz->inbuf, RZ_BUFFER_SIZE);
#endif
	ret = _read_gz_header(rz->inbuf, n, &ext_off, &ext_len);
	if(ret == 0){
		PLAIN_FILE:
		rz->in = n;
		rz->file_type = FILE_TYPE_PLAIN;
		memcpy(rz->outbuf, rz->inbuf, n);
		rz->buf_len = n;
		free(rz->stream);
		rz->stream = NULL;
		return rz;
	}
	rz->header_size = ret;
	ret = inflateInit2(rz->stream, -WINDOW_BITS);
	if(ret != Z_OK){ inflateEnd(rz->stream); goto PLAIN_FILE;}
	rz->stream->avail_in = n - rz->header_size;
	rz->stream->next_in  = rz->inbuf + rz->header_size;
	rz->stream->avail_out = RZ_BUFFER_SIZE;
	rz->stream->next_out  = rz->outbuf;
	rz->file_type = FILE_TYPE_GZ;
	rz->in = rz->header_size;
	rz->block_pos = rz->header_size;
	rz->next_block_pos = rz->header_size;
	rz->block_off = 0;
	if(ext_len < 7 || memcmp(rz->inbuf + ext_off, c, 4) != 0) return rz;
	if(((((unsigned char*)rz->inbuf)[ext_off + 5] << 8) | ((unsigned char*)rz->inbuf)[ext_off + 6]) != RZ_BLOCK_SIZE){
		fprintf(stderr, " -- WARNING: RZ_BLOCK_SIZE is not %d, treat source as gz file.  in %s -- %s:%d --\n", RZ_BLOCK_SIZE, __FUNCTION__, __FILE__, __LINE__);
		return rz;
	}
	rz->load_index = _load_index;
	rz->file_type = FILE_TYPE_RZ;
#ifdef _USE_KNETFILE
	if(knet_seek(fp, -16, SEEK_END) == -1){
#else
	if(lseek(fd, -16, SEEK_END) == -1){
#endif
		UNSEEKABLE:
		rz->seekable = 0;
		rz->index = NULL;
		rz->src_end = rz->end = 0x7FFFFFFFFFFFFFFFLL;
	} else {
		is_be = is_big_endian();
		rz->seekable = 1;
#ifdef _USE_KNETFILE
        knet_read(fp, &end, sizeof(int64_t));
#else
		read(fd, &end, sizeof(int64_t));
#endif        
		if(!is_be) rz->src_end = (int64_t)byte_swap_8((uint64_t)end);
		else rz->src_end = end;

#ifdef _USE_KNETFILE
		knet_read(fp, &end, sizeof(int64_t));
#else
		read(fd, &end, sizeof(int64_t));
#endif        
		if(!is_be) rz->end = (int64_t)byte_swap_8((uint64_t)end);
		else rz->end = end;
		if(n > rz->end){
			rz->stream->avail_in -= n - rz->end;
			n = rz->end;
		}
		if(rz->end > rz->src_end){
#ifdef _USE_KNETFILE
            knet_seek(fp, rz->in, SEEK_SET);
#else
			lseek(fd, rz->in, SEEK_SET);
#endif
			goto UNSEEKABLE;
		}
#ifdef _USE_KNETFILE
        knet_seek(fp, rz->end, SEEK_SET);
		if(knet_tell(fp) != rz->end){
			knet_seek(fp, rz->in, SEEK_SET);
#else
		if(lseek(fd, rz->end, SEEK_SET) != rz->end){
			lseek(fd, rz->in, SEEK_SET);
#endif
			goto UNSEEKABLE;
		}
#ifdef _USE_KNETFILE
		load_zindex(rz, fp);
		knet_seek(fp, n, SEEK_SET);
#else
		load_zindex(rz, fd);
		lseek(fd, n, SEEK_SET);
#endif
	}
	return rz;
}

#ifdef _USE_KNETFILE
RAZF* razf_dopen(int fd, const char *mode){
    if (strstr(mode, "r")) fprintf(stderr,"[razf_dopen] implement me\n");
    else if(strstr(mode, "w")) return razf_open_w(fd);
	return NULL;
}

RAZF* razf_dopen2(int fd, const char *mode)
{
    fprintf(stderr,"[razf_dopen2] implement me\n");
    return NULL;
}
#else
RAZF* razf_dopen(int fd, const char *mode){
	if(strstr(mode, "r")) return razf_open_r(fd, 1);
	else if(strstr(mode, "w")) return razf_open_w(fd);
	else return NULL;
}

RAZF* razf_dopen2(int fd, const char *mode)
{
	if(strstr(mode, "r")) return razf_open_r(fd, 0);
	else if(strstr(mode, "w")) return razf_open_w(fd);
	else return NULL;
}
#endif

static inline RAZF* _razf_open(const char *filename, const char *mode, int _load_index){
	int fd;
	RAZF *rz;
	if(strstr(mode, "r")){
#ifdef _USE_KNETFILE
        knetFile *fd = knet_open(filename, "r");
        if (fd == 0) {
            fprintf(stderr, "[_razf_open] fail to open %s\n", filename);
            return NULL;
        }
#else
#ifdef _WIN32
		fd = open(filename, O_RDONLY | O_BINARY);
#else
		fd = open(filename, O_RDONLY);
#endif
#endif
		if(fd < 0) return NULL;
		rz = razf_open_r(fd, _load_index);
	} else if(strstr(mode, "w")){
#ifdef _WIN32
		fd = open(filename, O_WRONLY | O_CREAT | O_TRUNC | O_BINARY, 0666);
#else
		fd = open(filename, O_WRONLY | O_CREAT | O_TRUNC, 0666);
#endif
		if(fd < 0) return NULL;
		rz = razf_open_w(fd);
	} else return NULL;
	return rz;
}

RAZF* razf_open(const char *filename, const char *mode){
	return _razf_open(filename, mode, 1);
}

RAZF* razf_open2(const char *filename, const char *mode){
	return _razf_open(filename, mode, 0);
}

int razf_get_data_size(RAZF *rz, int64_t *u_size, int64_t *c_size){
	int64_t n;
	if(rz->mode != 'r' && rz->mode != 'R') return 0;
	switch(rz->file_type){
		case FILE_TYPE_PLAIN:
			if(rz->end == 0x7fffffffffffffffLL){
#ifdef _USE_KNETFILE
				if(knet_seek(rz->x.fpr, 0, SEEK_CUR) == -1) return 0;
                n = knet_tell(rz->x.fpr);
				knet_seek(rz->x.fpr, 0, SEEK_END);
                rz->end = knet_tell(rz->x.fpr);
				knet_seek(rz->x.fpr, n, SEEK_SET);
#else
				if((n = lseek(rz->filedes, 0, SEEK_CUR)) == -1) return 0;
				rz->end = lseek(rz->filedes, 0, SEEK_END);
				lseek(rz->filedes, n, SEEK_SET);
#endif                
			}
			*u_size = *c_size = rz->end;
			return 1;
		case FILE_TYPE_GZ:
			return 0;
		case FILE_TYPE_RZ:
			if(rz->src_end == rz->end) return 0;
			*u_size = rz->src_end;
			*c_size = rz->end;
			return 1;
		default:
			return 0;
	}
}

static int _razf_read(RAZF* rz, void *data, int size){
	int ret, tin;
	if(rz->z_eof || rz->z_err) return 0;
	if (rz->file_type == FILE_TYPE_PLAIN) {
#ifdef _USE_KNETFILE
		ret = knet_read(rz->x.fpr, data, size);
#else
		ret = read(rz->filedes, data, size);
#endif        
		if (ret == 0) rz->z_eof = 1;
		return ret;
	}
	rz->stream->avail_out = size;
	rz->stream->next_out  = data;
	while(rz->stream->avail_out){
		if(rz->stream->avail_in == 0){
			if(rz->in >= rz->end){ rz->z_eof = 1; break; }
			if(rz->end - rz->in < RZ_BUFFER_SIZE){
#ifdef _USE_KNETFILE
				rz->stream->avail_in = knet_read(rz->x.fpr, rz->inbuf, rz->end -rz->in);
#else
				rz->stream->avail_in = read(rz->filedes, rz->inbuf, rz->end -rz->in);
#endif        
			} else {
#ifdef _USE_KNETFILE
				rz->stream->avail_in = knet_read(rz->x.fpr, rz->inbuf, RZ_BUFFER_SIZE);
#else
				rz->stream->avail_in = read(rz->filedes, rz->inbuf, RZ_BUFFER_SIZE);
#endif        
			}
			if(rz->stream->avail_in == 0){
				rz->z_eof = 1;
				break;
			}
			rz->stream->next_in = rz->inbuf;
		}
		tin = rz->stream->avail_in;
		ret = inflate(rz->stream, Z_BLOCK);
		rz->in += tin - rz->stream->avail_in;
		if(ret == Z_NEED_DICT || ret == Z_MEM_ERROR || ret == Z_DATA_ERROR){
			fprintf(stderr, "[_razf_read] inflate error: %d %s (at %s:%d)\n", ret, rz->stream->msg ? rz->stream->msg : "", __FILE__, __LINE__);
			rz->z_err = 1;
			break;
		}
		if(ret == Z_STREAM_END){
			rz->z_eof = 1;
			break;
		}
		if ((rz->stream->data_type&128) && !(rz->stream->data_type&64)){
			rz->buf_flush = 1;
			rz->next_block_pos = rz->in;
			break;
		}
	}
	return size - rz->stream->avail_out;
}

int razf_read(RAZF *rz, void *data, int size){
	int ori_size, i;
	ori_size = size;
	while(size > 0){
		if(rz->buf_len){
			if(size < rz->buf_len){
				for(i=0;i<size;i++) ((char*)data)[i] = ((char*)rz->outbuf + rz->buf_off)[i];
				rz->buf_off += size;
				rz->buf_len -= size;
				data += size;
				rz->block_off += size;
				size = 0;
				break;
			} else {
				for(i=0;i<rz->buf_len;i++) ((char*)data)[i] = ((char*)rz->outbuf + rz->buf_off)[i];
				data += rz->buf_len;
				size -= rz->buf_len;
				rz->block_off += rz->buf_len;
				rz->buf_off = 0;
				rz->buf_len = 0;
				if(rz->buf_flush){
					rz->block_pos = rz->next_block_pos;
					rz->block_off = 0;
					rz->buf_flush = 0;
				}
			}
		} else if(rz->buf_flush){
			rz->block_pos = rz->next_block_pos;
			rz->block_off = 0;
			rz->buf_flush = 0;
		}
		if(rz->buf_flush) continue;
		rz->buf_len = _razf_read(rz, rz->outbuf, RZ_BUFFER_SIZE);
		if(rz->z_eof && rz->buf_len == 0) break;
	}
	rz->out += ori_size - size;
	return ori_size - size;
}

int razf_skip(RAZF* rz, int size){
	int ori_size;
	ori_size = size;
	while(size > 0){
		if(rz->buf_len){
			if(size < rz->buf_len){
				rz->buf_off += size;
				rz->buf_len -= size;
				rz->block_off += size;
				size = 0;
				break;
			} else {
				size -= rz->buf_len;
				rz->buf_off = 0;
				rz->buf_len = 0;
				rz->block_off += rz->buf_len;
				if(rz->buf_flush){
					rz->block_pos = rz->next_block_pos;
					rz->block_off = 0;
					rz->buf_flush = 0;
				}
			}
		} else if(rz->buf_flush){
			rz->block_pos = rz->next_block_pos;
			rz->block_off = 0;
			rz->buf_flush = 0;
		}
		if(rz->buf_flush) continue;
		rz->buf_len = _razf_read(rz, rz->outbuf, RZ_BUFFER_SIZE);
		if(rz->z_eof || rz->z_err) break;
	}
	rz->out += ori_size - size;
	return ori_size - size;
}

static void _razf_reset_read(RAZF *rz, int64_t in, int64_t out){
#ifdef _USE_KNETFILE
	knet_seek(rz->x.fpr, in, SEEK_SET);
#else
	lseek(rz->filedes, in, SEEK_SET);
#endif
	rz->in  = in;
	rz->out = out;
	rz->block_pos = in;
	rz->next_block_pos = in;
	rz->block_off = 0;
	rz->buf_flush = 0;
	rz->z_eof = rz->z_err = 0;
	inflateReset(rz->stream);
	rz->stream->avail_in = 0;
	rz->buf_off = rz->buf_len = 0;
}

int64_t razf_jump(RAZF *rz, int64_t block_start, int block_offset){
	int64_t pos;
	rz->z_eof = 0;
	if(rz->file_type == FILE_TYPE_PLAIN){
		rz->buf_off = rz->buf_len = 0;
		pos = block_start + block_offset;
#ifdef _USE_KNETFILE
		knet_seek(rz->x.fpr, pos, SEEK_SET);
        pos = knet_tell(rz->x.fpr);
#else
		pos = lseek(rz->filedes, pos, SEEK_SET);
#endif
		rz->out = rz->in = pos;
		return pos;
	}
	if(block_start == rz->block_pos && block_offset >= rz->block_off) {
		block_offset -= rz->block_off;
		goto SKIP; // Needn't reset inflate
	}
	if(block_start  == 0) block_start = rz->header_size; // Automaticly revist wrong block_start
	_razf_reset_read(rz, block_start, 0);
	SKIP:
	if(block_offset) razf_skip(rz, block_offset);
	return rz->block_off;
}

int64_t razf_seek(RAZF* rz, int64_t pos, int where){
	int64_t idx;
	int64_t seek_pos, new_out;
	rz->z_eof = 0;
	if (where == SEEK_CUR) pos += rz->out;
	else if (where == SEEK_END) pos += rz->src_end;
	if(rz->file_type == FILE_TYPE_PLAIN){
#ifdef _USE_KNETFILE
		knet_seek(rz->x.fpr, pos, SEEK_SET);
        seek_pos = knet_tell(rz->x.fpr);
#else
		seek_pos = lseek(rz->filedes, pos, SEEK_SET);
#endif
		rz->buf_off = rz->buf_len = 0;
		rz->out = rz->in = seek_pos;
		return seek_pos;
	} else if(rz->file_type == FILE_TYPE_GZ){
		if(pos >= rz->out) goto SKIP;
		return rz->out;
	}
	if(pos == rz->out) return pos;
	if(pos > rz->src_end) return rz->out;
	if(!rz->seekable || !rz->load_index){
		if(pos >= rz->out) goto SKIP;
	}
	idx = pos / RZ_BLOCK_SIZE - 1;
	seek_pos = (idx < 0)? rz->header_size:(rz->index->cell_offsets[idx] + rz->index->bin_offsets[idx / RZ_BIN_SIZE]);
	new_out  = (idx + 1) * RZ_BLOCK_SIZE;
	if(pos > rz->out && new_out <= rz->out) goto SKIP;
	_razf_reset_read(rz, seek_pos, new_out);
	SKIP:
	razf_skip(rz, (int)(pos - rz->out));
	return rz->out;
}

uint64_t razf_tell2(RAZF *rz)
{
	/*
	if (rz->load_index) {
		int64_t idx, seek_pos;
		idx = rz->out / RZ_BLOCK_SIZE - 1;
		seek_pos = (idx < 0)? rz->header_size:(rz->index->cell_offsets[idx] + rz->index->bin_offsets[idx / RZ_BIN_SIZE]);
		if (seek_pos != rz->block_pos || rz->out%RZ_BLOCK_SIZE != rz->block_off)
			fprintf(stderr, "[razf_tell2] inconsistent block offset: (%lld, %lld) != (%lld, %lld)\n",
					(long long)seek_pos, (long long)rz->out%RZ_BLOCK_SIZE, (long long)rz->block_pos, (long long) rz->block_off);
	}
	*/
	return (uint64_t)rz->block_pos<<16 | (rz->block_off&0xffff);
}

int64_t razf_seek2(RAZF *rz, uint64_t voffset, int where)
{
	if (where != SEEK_SET) return -1;
	return razf_jump(rz, voffset>>16, voffset&0xffff);
}

void razf_close(RAZF *rz){
	if(rz->mode == 'w'){
#ifndef _RZ_READONLY
		razf_end_flush(rz);
		deflateEnd(rz->stream);
#ifdef _USE_KNETFILE
		save_zindex(rz, rz->x.fpw);
		if(is_big_endian()){
			write(rz->x.fpw, &rz->in, sizeof(int64_t));
			write(rz->x.fpw, &rz->out, sizeof(int64_t));
		} else {
			uint64_t v64 = byte_swap_8((uint64_t)rz->in);
			write(rz->x.fpw, &v64, sizeof(int64_t));
			v64 = byte_swap_8((uint64_t)rz->out);
			write(rz->x.fpw, &v64, sizeof(int64_t));
		}
#else
		save_zindex(rz, rz->filedes);
		if(is_big_endian()){
			write(rz->filedes, &rz->in, sizeof(int64_t));
			write(rz->filedes, &rz->out, sizeof(int64_t));
		} else {
			uint64_t v64 = byte_swap_8((uint64_t)rz->in);
			write(rz->filedes, &v64, sizeof(int64_t));
			v64 = byte_swap_8((uint64_t)rz->out);
			write(rz->filedes, &v64, sizeof(int64_t));
		}
#endif
#endif
	} else if(rz->mode == 'r'){
		if(rz->stream) inflateEnd(rz->stream);
	}
	if(rz->inbuf) free(rz->inbuf);
	if(rz->outbuf) free(rz->outbuf);
	if(rz->header){
		free(rz->header->extra);
		free(rz->header->name);
		free(rz->header->comment);
		free(rz->header);
	}
	if(rz->index){
		free(rz->index->bin_offsets);
		free(rz->index->cell_offsets);
		free(rz->index);
	}
	free(rz->stream);
#ifdef _USE_KNETFILE
    if (rz->mode == 'r')
        knet_close(rz->x.fpr);
    if (rz->mode == 'w')
        close(rz->x.fpw);
#else
	close(rz->filedes);
#endif
	free(rz);
}

#endif /* __RAZF_RJ_H */


/*source faidx.h and faidx.c*/
/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#ifndef FAIDX_H
#define FAIDX_H

/*!
  @header

  Index FASTA files and extract subsequence.

  @copyright The Wellcome Trust Sanger Institute.
 */

#include <ctype.h>

struct __faidx_t;
typedef struct __faidx_t faidx_t;


	/*!
	  @abstract   Build index for a FASTA or razip compressed FASTA file.
	  @param  fn  FASTA file name
	  @return     0 on success; or -1 on failure
	  @discussion File "fn.fai" will be generated.
	 */
	int fai_build(const char *fn);

	/*!
	  @abstract    Distroy a faidx_t struct.
	  @param  fai  Pointer to the struct to be destroyed
	 */
	void fai_destroy(faidx_t *fai);

	/*!
	  @abstract   Load index from "fn.fai".
	  @param  fn  File name of the FASTA file
	 */
	faidx_t *fai_load(const char *fn);

	/*!
	  @abstract    Fetch the sequence in a region.
	  @param  fai  Pointer to the faidx_t struct
	  @param  reg  Region in the format "chr2:20,000-30,000"
	  @param  len  Length of the region
	  @return      Pointer to the sequence; null on failure

	  @discussion The returned sequence is allocated by malloc family
	  and should be destroyed by end users by calling free() on it.
	 */
	char *fai_fetch(const faidx_t *fai, const char *reg, int *len);

	/*!
	  @abstract	   Fetch the number of sequences. 
	  @param  fai  Pointer to the faidx_t struct
	  @return	   The number of sequences
	 */
	int faidx_fetch_nseq(const faidx_t *fai);

	/*!
	  @abstract    Fetch the sequence in a region.
	  @param  fai  Pointer to the faidx_t struct
	  @param  c_name Region name
	  @param  p_beg_i  Beginning position number (zero-based)
	  @param  p_end_i  End position number (zero-based)
	  @param  len  Length of the region
	  @return      Pointer to the sequence; null on failure

	  @discussion The returned sequence is allocated by malloc family
	  and should be destroyed by end users by calling free() on it.
	 */
	char *faidx_fetch_seq(const faidx_t *fai, char *c_name, int p_beg_i, int p_end_i, int *len);

typedef struct {
	int32_t line_len, line_blen;
	int64_t len;
	uint64_t offset;
} faidx1_t;
KHASH_MAP_INIT_STR(s, faidx1_t)

struct __faidx_t {
	RAZF *rz;
        int n, m;
        int iter;
	char **name;
        char *filename;
	khash_t(s) *hash;
};

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

static inline void fai_insert_index(faidx_t *idx, const char *name, int len, int line_len, int line_blen, uint64_t offset)
{
	khint_t k;
	int ret;
	faidx1_t t;
	if (idx->n == idx->m) {
		idx->m = idx->m? idx->m<<1 : 16;
		idx->name = (char**)realloc(idx->name, sizeof(void*) * idx->m);
	}
	idx->name[idx->n] = strdup(name);
	k = kh_put(s, idx->hash, idx->name[idx->n], &ret);
	t.len = len; t.line_len = line_len; t.line_blen = line_blen; t.offset = offset;
	kh_value(idx->hash, k) = t;
	++idx->n;
}

faidx_t *fai_build_core(RAZF *rz)
{
	char c, *name;
	int l_name, m_name, ret;
	int line_len, line_blen, state;
	int l1, l2;
	faidx_t *idx;
	uint64_t offset;
	int64_t len;

	idx = (faidx_t*)calloc(1, sizeof(faidx_t));
	idx->hash = kh_init(s);
	name = 0; l_name = m_name = 0;
	len = line_len = line_blen = -1; state = 0; l1 = l2 = -1; offset = 0;
	while (razf_read(rz, &c, 1)) {
		if (c == '\n') { // an empty line
			if (state == 1) {
				offset = razf_tell(rz);
				continue;
			} else if ((state == 0 && len < 0) || state == 2) continue;
		}
		if (c == '>') { // fasta header
			if (len >= 0)
				fai_insert_index(idx, name, len, line_len, line_blen, offset);
			l_name = 0;
			while ((ret = razf_read(rz, &c, 1)) != 0 && !isspace(c)) {
				if (m_name < l_name + 2) {
					m_name = l_name + 2;
					kroundup32(m_name);
					name = (char*)realloc(name, m_name);
				}
				name[l_name++] = c;
			}
			name[l_name] = '\0';
			if (ret == 0) {
				fprintf(stderr, "[fai_build_core] the last entry has no sequence\n");
				free(name); fai_destroy(idx);
				return 0;
			}
			if (c != '\n') while (razf_read(rz, &c, 1) && c != '\n');
			state = 1; len = 0;
			offset = razf_tell(rz);
		} else {
			if (state == 3) {
				fprintf(stderr, "[fai_build_core] inlined empty line is not allowed in sequence '%s'.\n", name);
				free(name); fai_destroy(idx);
				return 0;
			}
			if (state == 2) state = 3;
			l1 = l2 = 0;
			do {
				++l1;
				if (isgraph(c)) ++l2;
			} while ((ret = razf_read(rz, &c, 1)) && c != '\n');
			if (state == 3 && l2) {
				fprintf(stderr, "[fai_build_core] different line length in sequence '%s'.\n", name);
				free(name); fai_destroy(idx);
				return 0;
			}
			++l1; len += l2;
			if (state == 1) line_len = l1, line_blen = l2, state = 0;
			else if (state == 0) {
				if (l1 != line_len || l2 != line_blen) state = 2;
			}
		}
	}
	fai_insert_index(idx, name, len, line_len, line_blen, offset);
	free(name);
	return idx;
}

void fai_save(const faidx_t *fai, FILE *fp)
{
	khint_t k;
	int i;
	for (i = 0; i < fai->n; ++i) {
		faidx1_t x;
		k = kh_get(s, fai->hash, fai->name[i]);
		x = kh_value(fai->hash, k);
#ifdef _WIN32
		fprintf(fp, "%s\t%d\t%ld\t%d\t%d\n", fai->name[i], (int)x.len, (long)x.offset, (int)x.line_blen, (int)x.line_len);
#else
		fprintf(fp, "%s\t%d\t%lld\t%d\t%d\n", fai->name[i], (int)x.len, (long long)x.offset, (int)x.line_blen, (int)x.line_len);
#endif
	}
}

faidx_t *fai_read(FILE *fp)
{
	faidx_t *fai;
	char *buf, *p;
	int len, line_len, line_blen;
#ifdef _WIN32
	long offset;
#else
	long long offset;
#endif
	fai = (faidx_t*)calloc(1, sizeof(faidx_t));
	fai->hash = kh_init(s);
	buf = (char*)calloc(0x10000, 1);
	while (!feof(fp) && fgets(buf, 0x10000, fp)) {
		for (p = buf; *p && isgraph(*p); ++p);
		*p = 0; ++p;
#ifdef _WIN32
		sscanf(p, "%d%ld%d%d", &len, &offset, &line_blen, &line_len);
#else
		sscanf(p, "%d%lld%d%d", &len, &offset, &line_blen, &line_len);
#endif
		fai_insert_index(fai, buf, len, line_len, line_blen, offset);
	}
	free(buf);
	return fai;
}

void fai_destroy(faidx_t *fai)
{
	int i;
	for (i = 0; i < fai->n; ++i) free(fai->name[i]);
	free(fai->name);
	free(fai->filename);
	kh_destroy(s, fai->hash);
	if (fai->rz) razf_close(fai->rz);
	free(fai);
}

int fai_build(const char *fn)
{
	char *str;
	RAZF *rz;
	FILE *fp;
	faidx_t *fai;
	str = (char*)calloc(strlen(fn) + 5, 1);
	sprintf(str, "%s.fai", fn);
	rz = razf_open(fn, "r");
	if (rz == 0) {
		fprintf(stderr, "[fai_build] fail to open the FASTA file %s\n",fn);
		free(str);
		return -1;
	}
	fai = fai_build_core(rz);
	razf_close(rz);
	fp = fopen(str, "wb");
	if (fp == 0) {
		fprintf(stderr, "[fai_build] fail to write FASTA index %s\n",str);
		fai_destroy(fai); free(str);
		return -1;
	}
	fai_save(fai, fp);
	fclose(fp);
	free(str);
	fai_destroy(fai);
	return 0;
}

#ifdef _USE_KNETFILE
FILE *download_and_open(const char *fn)
{
    const int buf_size = 1 * 1024 * 1024;
    uint8_t *buf;
    FILE *fp;
    knetFile *fp_remote;
    const char *url = fn;
    const char *p;
    int l = strlen(fn);
    for (p = fn + l - 1; p >= fn; --p)
        if (*p == '/') break;
    fn = p + 1;

    // First try to open a local copy
    fp = fopen(fn, "r");
    if (fp)
        return fp;

    // If failed, download from remote and open
    fp_remote = knet_open(url, "rb");
    if (fp_remote == 0) {
        fprintf(stderr, "[download_from_remote] fail to open remote file %s\n",url);
        return NULL;
    }
    if ((fp = fopen(fn, "wb")) == 0) {
        fprintf(stderr, "[download_from_remote] fail to create file in the working directory %s\n",fn);
        knet_close(fp_remote);
        return NULL;
    }
    buf = (uint8_t*)calloc(buf_size, 1);
    while ((l = knet_read(fp_remote, buf, buf_size)) != 0)
        fwrite(buf, 1, l, fp);
    free(buf);
    fclose(fp);
    knet_close(fp_remote);

    return fopen(fn, "r");
}
#endif

faidx_t *fai_load(const char *fn)
{
	char *str;
	FILE *fp;
	faidx_t *fai;
	str = (char*)calloc(strlen(fn) + 5, 1);
	sprintf(str, "%s.fai", fn);

#ifdef _USE_KNETFILE
    if (strstr(fn, "ftp://") == fn || strstr(fn, "http://") == fn)
    {
        fp = download_and_open(str);
        if ( !fp )
        {
            fprintf(stderr, "[fai_load] failed to open remote FASTA index %s\n", str);
            free(str);
            return 0;
        }
    }
    else
#endif
        fp = fopen(str, "rb");
	if (fp == 0) {
	        //fprintf(stderr, "[fai_load] build FASTA index.\n");
		fai_build(fn);
		fp = fopen(str, "rb");
		if (fp == 0) {
			fprintf(stderr, "[fai_load] fail to open FASTA index.\n");
			free(str);
			return 0;
		}
	}

	fai = fai_read(fp);
	fclose(fp);

	fai->rz = razf_open(fn, "rb");
	free(str);
	if (fai->rz == 0) {
		fprintf(stderr, "[fai_load] fail to open FASTA file.\n");
		return 0;
	}

	fai->filename = (char*)calloc(strlen(fn) + 1, 1);
	fai->iter = 0;
	strcpy(fai->filename, fn);

	return fai;
}

char *fai_fetch(const faidx_t *fai, const char *str, int *len)
{
	char *s, c;
	int i, l, k, name_end;
	khiter_t iter;
	faidx1_t val;
	khash_t(s) *h;
	int beg, end;

	beg = end = -1;
	h = fai->hash;
	name_end = l = strlen(str);
	s = (char*)malloc(l+1);
	// remove space
	for (i = k = 0; i < l; ++i)
		if (!isspace(str[i])) s[k++] = str[i];
	s[k] = 0; l = k;
	// determine the sequence name
	for (i = l - 1; i >= 0; --i) if (s[i] == ':') break; // look for colon from the end
	if (i >= 0) name_end = i;
	if (name_end < l) { // check if this is really the end
		int n_hyphen = 0;
		for (i = name_end + 1; i < l; ++i) {
			if (s[i] == '-') ++n_hyphen;
			else if (!isdigit(s[i]) && s[i] != ',') break;
		}
		if (i < l || n_hyphen > 1) name_end = l; // malformated region string; then take str as the name
		s[name_end] = 0;
		iter = kh_get(s, h, s);
		if (iter == kh_end(h)) { // cannot find the sequence name
			iter = kh_get(s, h, str); // try str as the name
			if (iter == kh_end(h)) {
				*len = 0;
			free(s); return 0;
			} else s[name_end] = ':', name_end = l;
		}
	} else iter = kh_get(s, h, str);
	if (iter == kh_end(h)) {
	  *len = 0;
	  free(s); return 0;
	}
	val = kh_value(h, iter);
	// parse the interval
	if (name_end < l) {
		for (i = k = name_end + 1; i < l; ++i)
			if (s[i] != ',') s[k++] = s[i];
		s[k] = 0;
		beg = atoi(s + name_end + 1);
		for (i = name_end + 1; i != k; ++i) if (s[i] == '-') break;
		end = i < k? atoi(s + i + 1) : val.len;
		if (beg > 0) --beg;
	} else beg = 0, end = val.len;
	if (beg >= val.len) beg = val.len;
	if (end >= val.len) end = val.len;
	if (beg > end) beg = end;
	free(s);

	// now retrieve the sequence
	l = 0;
	s = (char*)malloc(end - beg + 2);
	razf_seek(fai->rz, val.offset + beg / val.line_blen * val.line_len + beg % val.line_blen, SEEK_SET);
	while (razf_read(fai->rz, &c, 1) == 1 && l < end - beg && !fai->rz->z_err)
		if (isgraph(c)) s[l++] = c;
	s[l] = '\0';
	*len = l;
	return s;
}

int faidx_main(int argc, char *argv[])
{
	if (argc == 1) {
		fprintf(stderr, "Usage: faidx <in.fasta> [<reg> [...]]\n");
		return 1;
	} else {
		if (argc == 2) fai_build(argv[1]);
		else {
			int i, j, k, l;
			char *s;
			faidx_t *fai;
			fai = fai_load(argv[1]);
			if (fai == 0) return 1;
			for (i = 2; i != argc; ++i) {
				printf(">%s\n", argv[i]);
				s = fai_fetch(fai, argv[i], &l);
				for (j = 0; j < l; j += 60) {
					for (k = 0; k < 60 && k < l - j; ++k)
						putchar(s[j + k]);
					putchar('\n');
				}
				free(s);
			}
			fai_destroy(fai);
		}
	}
	return 0;
}

int faidx_fetch_nseq(const faidx_t *fai) 
{
	return fai->n;
}

char **faidx_all_keys(const faidx_t *fai)
{
  return fai->name;
}

int faidx_exists(const faidx_t *fai, char *c_name)
{
  khiter_t iter;

  iter = kh_get(s, fai->hash, c_name);
  if(iter == kh_end(fai->hash)){
    return 0;
  }
  else{
    return 1;
  }
}

int64_t faidx_seq_length(const faidx_t *fai, char *c_name)
{
  khiter_t iter;
  faidx1_t val;

  iter = kh_get(s, fai->hash, c_name);
  if(iter == kh_end(fai->hash)){
    return -1;
  }
  else{
    val = kh_value(fai->hash, iter);
    return val.len;
  }
}

char *faidx_fetch_seq(const faidx_t *fai, char *c_name, int p_beg_i, int p_end_i, int *len)
{
	int l;
	char c;
    khiter_t iter;
    faidx1_t val;
	char *seq=NULL;

    // Adjust position
    iter = kh_get(s, fai->hash, c_name);
    if(iter == kh_end(fai->hash)) return 0;
    val = kh_value(fai->hash, iter);
	if(p_end_i < p_beg_i) p_beg_i = p_end_i;
    if(p_beg_i < 0) p_beg_i = 0;
    else if(val.len <= p_beg_i) p_beg_i = val.len - 1;
    if(p_end_i < 0) p_end_i = 0;
    else if(val.len <= p_end_i) p_end_i = val.len - 1;

    // Now retrieve the sequence 
	l = 0;
	seq = (char*)malloc(p_end_i - p_beg_i + 2);
	razf_seek(fai->rz, val.offset + p_beg_i / val.line_blen * val.line_len + p_beg_i % val.line_blen, SEEK_SET);
	while (razf_read(fai->rz, &c, 1) == 1 && l < p_end_i - p_beg_i + 1)
		if (isgraph(c)) seq[l++] = c;
	seq[l] = '\0';
	*len = l;
	return seq;
}

char *faidx_fetch_header(const faidx_t *fai, char *c_name)
{
    int l;
    char c;
    khiter_t iter;
    faidx1_t val;
    char *head=NULL;
    char *tmp=NULL;
    int chunk=50;
    int size;
    int pos;

    /* Adjust position */
    iter = kh_get(s, fai->hash, c_name);
    if(iter == kh_end(fai->hash)) return 0;
    val = kh_value(fai->hash, iter);

    /* Now retrieve the header by walking backwards */
    l = 0;
    head = (char*)malloc(chunk);
    size = chunk;
    pos = val.offset - 1;
    
    while(pos > 0){
      /* Check if we need to grow string length */
      if (size <= l) {
	size += chunk;
	tmp = realloc(head, size);
	if (!tmp) {
	  free(head);
	  head = NULL;
	  break;
	}
	head = tmp;
      }

      /* Add character up to start of next line */
      pos--;
      razf_seek(fai->rz, pos, SEEK_SET);
      if(razf_read(fai->rz, &c, 1) == 1){
	if(c == '\n') break;
	head[l++] = c;
      }
      else{
	break;
      }
    }

    if(l == 0){
      free(head);
      return 0;
    }

    /* reorder the header (currently reversed) */
    head[l] = '\0';
    int i;
    for(i = 0; i < l; i++){
      int a = i;
      int b = l - 1 - i;
      if(b <= a) break;
      
      int buff = head[a];
      head[a] = head[b];
      head[b] = buff;
    }

    return head;
}

#endif /* FAIDX_H */


