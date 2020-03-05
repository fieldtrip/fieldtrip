# JAGA2FT

This application implements an interface from the Jaga16 device,
which streams data over WiFi using UDP, to the FieldTrip buffer.

The jaga16 will automatically connect to the WiFi network with SSID
jaganet and the password squeak!1

The UDB packets will be specifically broadcast to 192.168.8.100,
which is the IP address that should be granted by DHCP by the router,
or which should be configured as static address with router
192.168.8.1 and netmask 255.255.255.0.

