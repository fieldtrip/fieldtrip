function electrodesurface = cylinder_electrode(elec_coord,elec_normal,elec_rad,elec_thickness)

% CYLINDER_ELECTRODE returns a triangulatedsurface of a cylinder electrode

c0 = elec_coord + elec_normal; %make sure the bottom of the cylinder fully sticks into the surface it will be combined with afterwards
c1 = elec_coord - elec_thickness*elec_normal;

[electrodesurface.node,electrodesurface.face] = meshacylinder(c0,c1,elec_rad);
