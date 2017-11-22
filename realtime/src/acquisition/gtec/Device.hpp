//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#ifndef DEVICE_HPP_INCLUDED
#define DEVICE_HPP_INCLUDED

enum GDEVICE_TYPE 
{
    NO_TYPE, 
	GUSBAMP, 
	GHIAMP, 
	GNAUTILUS
};

struct Device 
{
    std::string name_;
    GDEVICE_TYPE type_;
    std::string full_type_;
    size_t state_;

    Device(std::string name, GDEVICE_TYPE type, size_t state) 
		: name_(name), type_(type), state_(state) 
	{
        if (type_ == GUSBAMP)
            full_type_ = "g.USBamp";
        else if (type_ == GHIAMP)
            full_type_ = "g.HIamp";
        else if (type_ == GNAUTILUS)
            full_type_ = "g.Nautilus";
    }

    Device()
		: name_(""), type_(NO_TYPE), full_type_(""), state_(0) 
	{
    }
};

#endif // DEVICE_HPP_INCLUDED
