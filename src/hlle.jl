#=
    Implements an adapted version of the hlle Riemann solver for the Shallow Water Equations
=#

import StaticArrays

# How to classify the two cells adjacent to the edge considered
@enum Cell_Pair_State begin
    UNDEFINED
    DRY__DRY
    WET__WET
    WET__DRY_INUNDATION
    WET__DRY_WALL
    WET__DRY_WALL_INUNDATION
    DRY__WET_INUNDATION
    DRY__WET_WALL
    DRY__WET_WALL_INUNDATION
end

# How to classify the Riemann problem on the edge considered
@enum RIEMANN_PROBLEM_TYPE begin
    DRY__SINGLE_RAREFACTION
    SINGLE_RAREFACTION__DRY
    SHOCK__SHOCK
    SHOCK__RAREFACTION
    RAREFACTION__SHOCK
    RAREFACTION__RAREFACTION
end

# The maximum water height a cell is allowed to have to still be considered dry
const DRY_TOL = 0.01

const ACCELERATION_GRAVITY = 9.81

const ZERO_TOL = 1.0e-7

# Can also be used to solve hv
function solve_riemann_hlle(
    h_left, h_right,
    hu_left, hu_right,
    b_left, b_right,
)
    # The quantities that have to be calculated in this routine
    flux_h_left = 0.0
    flux_h_right = 0.0
    flux_hu_left = 0.0
    flux_hu_right = 0.0
    max_wave_speed = 0.0

    # The speed of each cell (derived from the conserved quantities of momentum)
    u_left = 0.0
    u_right = 0.0

    ####
    # Determine dryness state
    ####
    wet_dry_state::Cell_Pair_State = UNDEFINED


    # compute speeds or set them to zero for the cell left of the edge
    if h_left >= DRY_TOL
        u_left = hu_left / h_left
    else
        b_left += h_left
        h_left = hu_left = u_left = 0.0
    end

    # compute speeds or set them to zero for the cell right of the edge
    if h_right >= DRY_TOL
        u_right = hu_right / h_right
    else
        b_right += h_right
        h_right = hu_right = u_right = 0.0
    end

    # determine wet/dry-state - try to start with the most frequent case
    if (h_left >= DRY_TOL)
        # the left cell is wet
        if h_right >= DRY_TOL
            # both cells are wet - expected to be the most frequent case
            wet_dry_state = WET__WET
        else
            # shore-line, left cell wet, right one is dry
            # check for inundation problems
            if h_left + b_left > b_right
                # dry cell (right) lies lower than wet cell (left)
                wet_dry_state = WET__DRY_INUNDATION
            else
                # dry cell (right) lies higher than wet cell (left)
                # -> remove middle state computation
                # -> assuming momentum is never high enough to overcome the
                # boundary on its own
                h_right = h_left
                u_right = - u_left
                hu_right = - hu_left
                b_right = b_left = 0.0
                wet_dry_state = WET__DRY_WALL
            end
        end
    else
        # the left cell is dry
        if h_right >= DRY_TOL
            # left cell dry, right cell wet
            # -> shoreline
            # -> check for inundation problems
            if h_right + b_right > b_left
                # dry cell (left) lies lower than wet cell (right)
                wet_dry_state = DRY__WET_INUNDATION
            else
                # dry cell (left) lies higher than wet cell (right)
                # -> remove middle state computation
                # -> assuming momentum is never high enough to overcome the
                # boundary on its own
                h_left = h_right
                u_left = - u_right
                hu_left = - hu_right
                b_left = b_right = 0.0
                wet_dry_state = DRY__WET_WALL
            end
        end
    end
    
    # compute some commonly used term to avoid expensive evaluations of the sqrt
    sqrt_h_left = sqrt(h_left)
    sqrt_h_right = sqrt(h_right)


    # compute eigenvalues of the jacobian matrix (== characteristic speeds)
    # Actually there are three since we have 3 conserved quantities but one
    # speed will always be 0
    characteristic_speeds = StaticArrays.SA_F64[
        u_left - sqrt(ACCELERATION_GRAVITY * h_left),
        u_right + sqrt(ACCELERATION_GRAVITY * h_right),
    ]

    # compute the "Roe" averaged quantities
    h_roe = 0.5 * (h_right + h_left)
    u_roe = (u_left * sqrt_h_left + u_right * sqrt_h_right) /
        (sqrt_h_left + sqrt_h_right)
    
    # compute the "Roe" speeds (i.e. quantities similar to the characteristic
    # speeds)
    roe_speeds = StaticArrays.SA_F64[
        u_roe - sqrt(ACCELERATION_GRAVITY * h_roe),
        u_roe + sqrt(ACCELERATION_GRAVITY * h_roe),
    ]

    # compute the einfeldt speeds based on the situation
    extremal_einfeldt_speeds = StaticArrays.MVector{2, Float64}(undef)
    if wet_dry_state == WET__WET || wet_dry_state == WET__DRY_WALL || wet_dry_state == DRY__WET_WALL
        extremal_einfeldt_speeds[1] =
            min(characteristic_speeds[1], roe_speeds[1])
        extremal_einfeldt_speeds[2] =
            max(characteristic_speeds[2], roe_speeds[2])
    elseif h_left < DRY_TOL
        # contains cases DRY__WET_INUNDATION and DRY__WET_WALL_INUNDATION
        extremal_einfeldt_speeds[1] = roe_speeds[1]
        extremal_einfeldt_speeds[2] =
            max(characteristic_speeds[2], roe_speeds[2])
    elseif h_right < DRY_TOL
        # contains cases WET__DRY_INUNDATION and WET__DRY_WALL_INUNDATION
        extremal_einfeldt_speeds[1] =
            min(characteristic_speeds[1], roe_speeds[1])
        extremal_einfeldt_speeds[2] = roe_speeds[2]
    end

    # compute the middle state according to the hll algorithm
    hll_height_middle = max(
        (hu_left - hu_right +
        extremal_einfeldt_speeds[2] * h_right -
        extremal_einfeldt_speeds[1] * h_left) /
        (extremal_einfeldt_speeds[2] - extremal_einfeldt_speeds[1]),
        0.0
    )

    # define the approximate eigenvalues of the Riemann problem
    eigenvalues = StaticArrays.SA_F64[
        extremal_einfeldt_speeds[1],
        0.5 * (extremal_einfeldt_speeds[1] + extremal_einfeldt_speeds[2]),
        extremal_einfeldt_speeds[2],
    ]

    # define the approximate eigenvectors of the Riemann problem
    # ugly construct with transpose necessary since julia is column-major
    eigenvectors = StaticArrays.SA_F64[
        1.0 0.0 1.0;
        eigenvalues[1] 0.0 eigenvalues[3];
        eigenvalues[1]^2 1.0 eigenvalues[3]^2 
    ]
    #=
        StaticArrays.SA_F64[
            1.0,
            0.0,
            1.0,
        ],
        StaticArrays.SA_F64[
            eigenvalues[1],
            0.0,
            eigenvalues[3],
        ],
        StaticArrays.SA_F64[
            eigenvalues[1] * eigenvalues[1],
            1.0,
            eigenvalues[3] * eigenvalues[3],
        ],
    ]
    =#

    # compute the jump in state, this will be the rhs of system of Equations
    right_hand_side = StaticArrays.MVector(
        h_right - h_left,
        hu_right - hu_left,
        (hu_right * u_right + 0.5 * ACCELERATION_GRAVITY * h_right^2) -
            (hu_left * u_left + 0.5 * ACCELERATION_GRAVITY * h_left^2)
    )

    # compute the steady-state wave
    steady_state_wave = StaticArrays.MVector(
        - (b_right - b_left),
        - 0.5 * ACCELERATION_GRAVITY * (h_left + h_right) * (b_right - b_left),
    )

    # preserve the depth-positivity i.e. find out whether the wave speed
    # directions is subsonic (the two extremal waves both have their own
    # quadrant) or supersonic (the two extremal waves both lie in the same
    # quadrant either both forward-facing or backward-facing)
    if eigenvalues[1] < - ZERO_TOL && eigenvalues[3] > ZERO_TOL
        # subsonic (spreading wave pattern)
        steady_state_wave[1] = max(
            steady_state_wave[1],
            hll_height_middle * (eigenvalues[3] - eigenvalues[1]) /
                eigenvalues[1]
        )
        steady_state_wave[1] = min(
            steady_state_wave[1],
            hll_height_middle * (eigenvalues[3] - eigenvalues[1]) /
                eigenvalues[3]
        )
    elseif eigenvalues[1] > ZERO_TOL
        # supersonic to the right
        steady_state_wave[1] = max(
            steady_state_wave[1],
            - h_left
        )
        steady_state_wave[1] = min(
            steady_state_wave[1],
            hll_height_middle * (eigenvalues[3] - eigenvalues[1]) /
                eigenvalues[1]
        )
    elseif eigenvalues[3] < - ZERO_TOL
        # supersonic to the left
        steady_state_wave[1] = max(
            steady_state_wave[1],
            hll_height_middle * (eigenvalues[3] - eigenvalues[1]) /
                eigenvalues[3]
        )
        steady_state_wave[1] = min(
            steady_state_wave[1],
            h_right
        )
    end

    # limit the effect of the source term
    steady_state_wave[2] = min(
        steady_state_wave[2],
        ACCELERATION_GRAVITY * max(
            - h_left * (b_right - b_left),
            - h_right * (b_right - b_left)
        )
    )
    steady_state_wave[2] = max(
        steady_state_wave[2],
        ACCELERATION_GRAVITY * min(
            - h_left * (b_right - b_left),
            - h_right * (b_right - b_left)
        )
    )

    # Apply effects of the steady-state wave
    right_hand_side[1] -= steady_state_wave[1]
    # right_hand_side[2]: no source term
    right_hand_side[3] -= steady_state_wave[2]

    # solve the linear system of equations 
    # eigenvalue_matrix * beta = right_hand_side 
    # with the unkowns beta which are the later coefficients to
    # calculate the f-waves
    #
    # the solution process is basically based on gaussian elimination
    beta = StaticArrays.MVector{3, Float64}(undef)
    beta[1] = (eigenvalues[3] * right_hand_side[1] - right_hand_side[2]) /
        (eigenvalues[3] - eigenvalues[1])
    beta[3] = (- eigenvalues[1] * right_hand_side[1] + right_hand_side[2]) /
        (eigenvalues[3] - eigenvalues[1])
    beta[2] = right_hand_side[3] -
        eigenvalues[1] * eigenvalues[1] * beta[1] -
        eigenvalues[3] * eigenvalues[3] * beta[3]
    
    

    # compute f-waves and wave speeds
    f_waves = StaticArrays.MMatrix{3, 2, Float64}(undef)
    f_waves .= 0.0
    wave_speeds = StaticArrays.MVector{3, Float64}(undef)
    wave_speeds .= 0.0

    if wet_dry_state == WET__DRY_WALL
        # zero ghost updates (wall boundary)
        # care about the left-going wave (index 1) only
        f_waves[1, 1] = beta[1] * eigenvectors[2, 1]
        f_waves[1, 2] = beta[1] * eigenvectors[3, 1]

        wave_speeds[1] = eigenvalues[1]
    elseif wet_dry_state == DRY__WET_WALL
        # zero ghost updates (wall boundary)
        # care about the right-going wave (index 3) only
        f_waves[3, 1] = beta[3] * eigenvectors[2, 3]
        f_waves[3, 2] = beta[3] * eigenvectors[3, 3]

        wave_speeds[3] = eigenvalues[3]
    else
        # default computation
        for wave_number in 1:3
            f_waves[wave_number, 1] = beta[wave_number] *
                eigenvectors[2, wave_number]
            f_waves[wave_number, 2] = beta[wave_number] * 
                eigenvectors[3, wave_number]
        end
        # TODO check if deepcopy might be necessary
        wave_speeds = eigenvalues
    end

    # Use the f-waves to calculate the numerical fluxes into each of the two
    # adjacent cells for each of the occuring wave
    for wave_number in 1:3
        if wave_speeds[wave_number] < - ZERO_TOL
            # left going
            flux_h_left += f_waves[wave_number, 1]
            flux_hu_left += f_waves[wave_number, 2]
        elseif wave_speeds[wave_number] > ZERO_TOL
            # right going
            flux_h_right += f_waves[wave_number, 1]
            flux_hu_right += f_waves[wave_number, 2]
        else
            # Case shouldn't happen mathematically
            # but does
            #println("Entered case that shouldn't happen")
            flux_h_left += 0.5 * f_waves[wave_number, 1]
            flux_hu_left += 0.5 * f_waves[wave_number, 2]

            flux_h_right += 0.5 * f_waves[wave_number, 1]
            flux_hu_right += 0.5 * f_waves[wave_number, 2]
        end
    end

    max_wave_speed = max(abs(wave_speeds[1]), abs(wave_speeds[2]), wave_speeds[3])

    return flux_h_left, flux_h_right, flux_hu_left, flux_hu_right, max_wave_speed
end