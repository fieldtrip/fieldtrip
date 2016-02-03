
%% Version 2016-16-01: Major update
% *New network models*
%
% * generative_model.m: Implements more than 10 generative network models.
% * evaluate_generative_model.m: Implements and evaluates the accuracy of
% more than 10 generative network models.
% * demo_generative_models_geometric.m and
% demo_generative_models_neighbors.m: Demonstrate the capabilities of the 
% new generative model functions.
%
% *New network measures*
%
% * clustering_coef_wu_sign.m: Multiple generalizations of the clustering
% coefficient for networks with positive and negative weights.
% * core_periphery_dir.m: Optimal core structure and core-ness statistic.
% * gateway_coef_sign.m: Gateway coefficient (a variant of the
% participation coefficient) for networks with positive and negative
% weights.
% * local_assortativity_sign.m: Local (nodal) assortativity for networks
% with positive and negative weights.
% * randmio_dir_signed.m: Random directed graph with preserved signed in-
% and out- degree distribution.
%
% *Removed network measures*
%
% * modularity_louvain_und_sign.m, modularity_finetune_und_sign.m: This
% functionality is now provided by community_louvain.m.
% * modularity_probtune_und_sign.m: Similar functionality is provided by 
% consensus_und.m
%
% *Bug fixes and/or code improvements and/or documentation improvements*
% 
% * charpath.m: Changed default behavior, such that infinitely long paths
% (i.e. paths between disconnected nodes) are now included in computations
% by default, but may be excluded manually.
% * community_louvain.m: Included generalization for negative weights, 
% enforced binary network input for Potts-model Hamiltonian, streamlined
% code.
% * eigenvector_centrality_und.m: Ensured the use of leading eigenvector
% for computations of eigenvector centrality.
% * modularity_und.m, modularity_dir.m: Enforced single node moves during
% fine-tuning step.
% * null_model_und_sign.m and null_model_dir_sign.m: Fixed preservation
% of negative degrees in sparse networks with negative weights.
% * randmio_und_signed.m: Now allows unbiased exploration of all network
% configurations.
% * transitivity_bd.m, transitivity_wu.m, transitivity_wd.m: removed tests
% for absence of nodewise 3-cycles. Expanded documentation.
% * clustering_coef_wu.m, clustering_coef_wd.m: Expanded documentation.
% * motif3-m and motif4-m functions: Expanded documentation.
% * rich_club_wu.m, rich_club_wd.m. Expanded documentation.
% 
% *Cosmetic and MATLAB code analyzer (mlint) improvements to many other functions*
%
%% Version 2015-25-01: Major update
% Includes two new community-detection scripts and multiple improvements 
% 
% * New community detection scripts: 1. community_louvain.m (supersedes
% modularity_louvain.m and modularity_finetune.m scripts); 2.
% link_communities.m.
% * added autofix flag to weight_conversion.m for fixing common weight
% problems.
% * other function improvements: participation_coef.m, charpath.m,
% reorder_mod.m.
% * bug fixes: modularity_finetune_und_sign.m,
% modularity_probtune_und_sign.m, threshold_proportional.m 
% * changed help files: assortativity_wei.m, distance_wei.m
% 
% 
%% Version 2014-04-05: Minor update
%
% * consensus_und.m is now a self-contained function
% * headers in charpath.m and in threshold_proportional.m have been corrected
