// Stub DBReader.h - minimal declaration to satisfy includes from extracted MMseqs2 code.
// Only provides the static utility methods actually used (getExtendedDbtype, setExtendedDbtype).
#ifndef DBREADER_H
#define DBREADER_H

#include <cstdint>
#include "Parameters.h"

template <typename T>
class DBReader {
public:
    static inline uint16_t getExtendedDbtype(int dbtype) {
        return (uint16_t)((uint32_t)dbtype >> 16) & 0x7FFE;
    }

    static inline int setExtendedDbtype(int dbtype, uint16_t extended) {
        return dbtype | ((extended & 0x7FFE) << 16);
    }
};

#endif // DBREADER_H
