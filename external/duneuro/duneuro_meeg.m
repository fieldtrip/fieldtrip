classdef duneuro_meeg < handle
    properties (Hidden = true)
        cpp_handle;
        constructor_arguments;
        source_model;
        electrodes;
        coils;
        projections;
    end
    methods
        % Constructor
        function this = duneuro_meeg(config)
            this.cpp_handle = duneuro_matlab('create', config);
            this.constructor_arguments = config;
            this.source_model = [];
            this.electrodes = [];
            this.coils = [];
            this.projections = [];
        end
        % Destructor
        function delete(this)
            duneuro_matlab('delete', this.cpp_handle);
        end
        function solve_eeg_forward(this, dipole, func, config)
            duneuro_matlab('solve_eeg_forward', this.cpp_handle, dipole, func.cpp_handle, config);
        end
        function solution = solve_meg_forward(this, func, config)
            solution = duneuro_matlab('solve_meg_forward', this.cpp_handle, func.cpp_handle, config);
        end
        function matrix = compute_eeg_transfer_matrix(this, config)
            matrix = duneuro_matlab('compute_eeg_transfer_matrix', this.cpp_handle, config);
        end
        function matrix = compute_meg_transfer_matrix(this, config)
            matrix = duneuro_matlab('compute_meg_transfer_matrix', this.cpp_handle, config);
        end
        function set_electrodes(this, electrodes, config)
            duneuro_matlab('set_electrodes', this.cpp_handle, electrodes, config);
            this.electrodes.electrodes = electrodes;
            this.electrodes.config = config;
        end
        function electrodes = get_projected_electrodes(this)
            electrodes = duneuro_matlab('get_projected_electrodes', this.cpp_handle);
        end
        function set_coils_and_projections(this, coils, projections)
            duneuro_matlab('set_coils_and_projections', this.cpp_handle, coils, projections);
            this.coils = coils;
            this.projections = projections;
        end
        function solution = evaluate_at_electrodes(this, func)
            solution = duneuro_matlab('evaluate_at_electrodes', this.cpp_handle, func.cpp_handle);
        end
        function solution = apply_eeg_transfer(this, transfer_matrix, dipoles, config)
            solution = duneuro_matlab('apply_eeg_transfer', this.cpp_handle, transfer_matrix, dipoles, config);
        end
        function solution = apply_meg_transfer(this, transfer_matrix, dipoles, config)
            solution = duneuro_matlab('apply_meg_transfer', this.cpp_handle, transfer_matrix, dipoles, config);
        end
        function write(this, func, config)
            duneuro_matlab('write', this.cpp_handle, func.cpp_handle, config);
        end
        function write_mesh(this, config)
            duneuro_matlab('write', this.cpp_handle, config);
        end
        function print_citations(this)
            duneuro_matlab('print_citations', this.cpp_handle);
        end
        function s = saveobj(this)
            s.constructor_arguments = this.constructor_arguments;
            s.source_model = this.source_model;
            s.electrodes = this.electrodes;
            s.coils = this.coils;
            s.projections = this.projections;
        end
    end
    methods(Static)
        function obj = loadobj(s)
            obj = duneuro_meeg(s.constructor_arguments);
            if ~isempty(s.source_model)
                obj.set_source_model(s.source_model);
            end
            if ~isempty(s.electrodes)
                obj.set_electrodes(s.electrodes.electrodes, s.electrodes.config);
            end
            if ~isempty(s.coils) && ~isempty(s.projections)
                obj.set_coils_and_projections(s.coils, s.projections);
            end
        end
    end
end
