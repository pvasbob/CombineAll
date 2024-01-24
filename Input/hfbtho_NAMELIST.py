import math


hfbtho_NAMELIST = {
    '&hfbtho_general': {
        'number_of_shells': 6,
        'oscillator_length': -1.0,
        'basis_deformation': +0.15,
        'proton_number': 12,
        'neutron_number': 12,
        'type_of_calculation': 1
    },
    '&hfbtho_initial': {
        'beta2_deformation': +0.15,
        'beta3_deformation': 0.0,
        'beta4_deformation': 0.0
    },
    '&hfbtho_iterations': {
        'number_iterations': 300,
        'accuracy': 1e-06,
        'restart_file': 1
    },
    '&hfbtho_functional': {
        'functional': 'SKM*',
        'add_initial_pairing': False,
        'type_of_coulomb': 2
    },
    '&hfbtho_pairing': {
        'user_pairing': False,
        'vpair_n': -250.0,
        'vpair_p': -250.0,
        'pairing_cutoff': 60.0,
        'pairing_feature': 0.5
    },
    '&hfbtho_constraints': {
        'lambda_values': [1, 2, 3, 4, 5, 6, 7, 8],
        'lambda_active': [0, -1, 0, 0, 0, 0, 0, 0],
        'expectation_values': [0.0, 0.5927836524190626, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    },
    '&hfbtho_blocking': {
        'proton_blocking': [0, 0, 0, 0, 0],
        'neutron_blocking': [0, 0, 0, 0, 0]
    },
    '&hfbtho_projection': {
        'switch_to_tho': 0,
        'projection_is_on': 0,
        'gauge_points': 1,
        'delta_z': 0,
        'delta_n': 0
    },
    '&hfbtho_temperature': {
        'set_temperature': False,
        'temperature': 0.0
    },
    '&hfbtho_features': {
        'collective_inertia': False,
        'fission_fragments': False,
        'pairing_regularization': False,
        'automatic_basis': False,
        'localization_functions': False
    },
    '&hfbtho_tddft': {
        'filter': False,
        'fragment_properties': False,
        'real_z': 49.7,
        'real_n': 71.2
    },
    '&hfbtho_neck': {
        'set_neck_constrain': False,
        'neck_value': 13.0
    },
    '&hfbtho_debug': {
        'number_gauss': 20,
        'number_laguerre': 20,
        'number_legendre': 40,
        'compatibility_hfodd': False,
        'number_states': 500,
        'force_parity': True,
        'print_time': 0
    },
    '&hfbtho_restoration': {
        'pnp_is_on': 0,
        'number_of_gauge_points': 9,
        'delta_neutrons': 6,
        'delta_protons': 6,
        'amp_is_on': 0,
        'number_of_rotational_angles': 27,
        'maximal_angular_momentum': 10
    }
}

# Updat constraints so that it satisfiy:
A =  hfbtho_NAMELIST['&hfbtho_general']['proton_number'] + hfbtho_NAMELIST['&hfbtho_general']['neutron_number']
hfbtho_NAMELIST['&hfbtho_constraints']['expectation_values'][1] = (math.sqrt(5 / math.pi) / 100) * hfbtho_NAMELIST['&hfbtho_general']['basis_deformation'] * A**(5/3)



