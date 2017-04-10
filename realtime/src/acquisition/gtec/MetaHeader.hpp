//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#ifndef METAHEADER_HPP_INCLUDED
#define METAHEADER_HPP_INCLUDED

#include <stdint.h>

static const uint64_t GDSMSG = 0x3E47534D5344473C; // <GDSMSG> in reverse byte order
static uint64_t current_package_id = 0; // consecutive packge number

#pragma pack(1)
#ifdef __cplusplus
extern "C" {
#endif

struct MetaHeader {
    uint64_t header_identifier_;
    uint64_t package_id_;
    uint8_t acknowledge_;
    uint64_t size_;
    uint8_t crc_enabled_;
    uint16_t crc_value_;

    explicit MetaHeader(uint64_t size) 
		: header_identifier_(GDSMSG), 
		package_id_(current_package_id++), 
		acknowledge_(0), 
		size_(size), 
		crc_enabled_(false), 
		crc_value_(0) 
	{
    }

    MetaHeader() 
		: header_identifier_(GDSMSG), 
		package_id_(current_package_id++), 
		acknowledge_(0), 
		size_(0), 
		crc_enabled_(false), 
		crc_value_(0) 
	{
    }
};

#ifdef __cplusplus
}
#endif
#pragma pack()

#endif // METAHEADER_HPP_INCLUDED
