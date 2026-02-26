#ifndef __UTILS__
#define __UTILS__

#define CACHE_BLOCK_SIZE 0

#define NUMBER_CORES 0

//A -> 65 -> 0b1000001 -> 0x41
//C -> 67 -> 0b1000011 -> 0x43
//G -> 71 -> 0b1000111 -> 0x47
//T -> 84 -> 0b1010100 -> 0x54
//N -> 78 -> 0b1001110 -> 0x4E

//a -> 97  -> 0b1100001 -> 0x61
//c -> 99  -> 0b1100011 -> 0x63
//g -> 103 -> 0b1100111 -> 0x67
//t -> 116 -> 0b1101110 -> 0x74
//n -> 110 -> 0b1110100 -> 0x6E

//char r;
//int gc_count = 0;
//gc_count += ((r & 0x3) == 0x3);

#define gc_matching(char_read_) \
    ((char_read_ & 0x3) == 0x3)

#endif