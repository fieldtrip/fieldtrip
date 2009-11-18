% BayesBrain is a Matlab graphical model toolbox. It is specifically tailored
% to work together with FieldTrip (R) for the purpose of single-trial analysis
% in cognitive neuroscience but can be used as a standalone toolbox as well
%
% Copyrights (C) 2008
% Intelligent Systems (http://www.ru.is/ml)
% FC Donders Centre (http://www.ru.nl/fcdonders/)
% Radboud University Nijmegen, The Netherlands
%
% Most functions in this toolbox are licensed under the GNU General
% Public license (GPL). Unauthorised copying and distribution of functions 
% that are not explicitely covered by the GPL is not allowed!
%
% The functions in this toolbox are copyrighted by their authors:
%   Marcel van Gerven, Intelligent Systems, FCDC
%
% The goal of BayesBrain is to provide a toolbox that allows learning of and
% inference with graphical models in cognitive neuroscience. The toolbox
% can be used to construct static and dynamic models.
%
% Some of the functionality relies on external toolboxes. E.g., inference
% with hugin_ie requires the commercial Hugin library, which is not
% included in the distribution.
%
% SEE ALSO
% graphicalmodel.m
% inference_engine.m
% cpd.m
% pot.m
% graph.m
%
% TO DO
% - construct other inference algorithms
% - make other inference algorithms suitable for DBNs (e.g., exhaustive i.e.)
% - create smoothing engine
% - make code more transparent
% - can we remove either cliques or potdom in jtree object?
% - Bug report Arjen Hommersom
