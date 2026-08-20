# Miranda & Rafikov (2019/2020) locally-cooled linear disc-mode MASTER-EQUATION COEFFICIENTS.
#
# Ported verbatim from the legacy `DiscModes/master_equation_cool_newt.jl`, itself the machine
# translation of the authors' `planet-disk-modes-cooling/master_equation_cool.py` (commit 948d41c)
# — MR20 Eqs. (9), (24)-(26). Pure arithmetic (only `**` → `^`), TYPE-GENERIC: the enthalpy master
# equation δh″ = Ψ − C₁δh′ − C₀δh split into real and imaginary parts.
#
# This is GENERATED code: it is transcribed, never hand-edited. The gate is
# `test/Modes/test_newtonian_planet.jl` (baked one-shot capture of the legacy package).
function master_equation_cool(r, gamma, cs, Omega, dOmegadr, Sigma, dSigmadr, LTinv, dLTinvdr, LSinv, dLSinvdr, omega_tilde, tc, dtcdr, D, dDdr, Nr2, dNr2dr, pole_disp, Phi, dPhidr, d2Phidr2, m)

    C1R = dSigmadr / Sigma + (
        D ^ 2
        * (
            1
            + 2 * gamma * pole_disp * tc
            + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
        )
        * (
            1
            + 2 * gamma * pole_disp * tc
            + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
            - LTinv * r * (1 + gamma * pole_disp * tc)
        )
        + gamma
        * Nr2
        * tc
        * (
            -(
                r
                * (
                    dDdr * pole_disp
                    + (dDdr + dNr2dr) * gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc
                )
                * (
                    1
                    + 2 * gamma * pole_disp * tc
                    + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                )
            )
            + gamma
            * Nr2
            * (
                -(dtcdr * (omega_tilde ^ 2 + pole_disp ^ 2) * r)
                + (
                    pole_disp ^ 2
                    + dOmegadr * m * omega_tilde * r
                    - dtcdr * gamma * pole_disp ^ 3 * r
                    - LTinv * (omega_tilde ^ 2 + pole_disp ^ 2) * r
                    + omega_tilde ^ 2 * (1 - dtcdr * gamma * pole_disp * r)
                )
                * tc
                + gamma
                * pole_disp
                * (
                    -(LTinv * (omega_tilde ^ 2 + pole_disp ^ 2) * r)
                    + 2
                    * (
                        omega_tilde ^ 2
                        + pole_disp ^ 2
                        + dOmegadr * m * omega_tilde * r
                    )
                )
                * tc ^ 2
                + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) ^ 2 * tc ^ 3
            )
        )
        + D
        * (
            -(
                r
                * (
                    1
                    + 2 * gamma * pole_disp * tc
                    + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                )
                * (
                    dDdr
                    + (2 * dDdr + dNr2dr) * gamma * pole_disp * tc
                    + (dDdr + dNr2dr)
                    * gamma ^ 2
                    * (omega_tilde ^ 2 + pole_disp ^ 2)
                    * tc ^ 2
                )
            )
            + gamma
            * Nr2
            * (
                -(dtcdr * pole_disp * r)
                - 2
                * (
                    -pole_disp
                    + LTinv * pole_disp * r
                    + dtcdr * gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * r
                )
                * tc
                + gamma
                * (
                    2 * dOmegadr * m * omega_tilde * r
                    - 2 * LTinv * (omega_tilde ^ 2 + 2 * pole_disp ^ 2) * r
                    + omega_tilde ^ 2 * (2 - dtcdr * gamma * pole_disp * r)
                    + pole_disp ^ 2 * (6 - dtcdr * gamma * pole_disp * r)
                )
                * tc ^ 2
                + 2
                * gamma ^ 2
                * pole_disp
                * (
                    3 * omega_tilde ^ 2
                    + 3 * pole_disp ^ 2
                    + dOmegadr * m * omega_tilde * r
                    - LTinv * (omega_tilde ^ 2 + pole_disp ^ 2) * r
                )
                * tc ^ 3
                + 2 * gamma ^ 3 * (omega_tilde ^ 2 + pole_disp ^ 2) ^ 2 * tc ^ 4
            )
        )
    ) / (
        r
        * (
            1
            + 2 * gamma * pole_disp * tc
            + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
        )
        * (
            gamma ^ 2 * Nr2 ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
            + 2
            * D
            * gamma
            * Nr2
            * tc
            * (pole_disp + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc)
            + D ^ 2
            * (
                1
                + 2 * gamma * pole_disp * tc
                + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
            )
        )
    )

    C1I = (
        gamma
        * (
            -(LTinv * omega_tilde * tc)
            + (
                dOmegadr * gamma ^ 2 * m * Nr2 * (D + Nr2) * omega_tilde ^ 2 * tc ^ 3
                + gamma ^ 2
                * omega_tilde ^ 3
                * tc ^ 2
                * (
                    -(Nr2 * (dtcdr * Nr2 + dDdr * tc))
                    + D * (-(dtcdr * Nr2) + dNr2dr * tc)
                )
                - dOmegadr
                * m
                * Nr2
                * tc
                * (1 + gamma * pole_disp * tc)
                * (gamma * Nr2 * pole_disp * tc + D * (1 + gamma * pole_disp * tc))
                + omega_tilde
                * (
                    D
                    * (1 + gamma * pole_disp * tc)
                    * (
                        dtcdr * Nr2
                        + (dNr2dr - dtcdr * gamma * Nr2 * pole_disp) * tc
                        + dNr2dr * gamma * pole_disp * tc ^ 2
                    )
                    - Nr2
                    * tc
                    * (
                        dtcdr * gamma ^ 2 * Nr2 * pole_disp ^ 2 * tc
                        + dDdr * (1 + gamma * pole_disp * tc) ^ 2
                    )
                )
            )
            / (
                gamma ^ 2 * Nr2 ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                + 2
                * D
                * gamma
                * Nr2
                * tc
                * (pole_disp + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc)
                + D ^ 2
                * (
                    1
                    + 2 * gamma * pole_disp * tc
                    + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                )
            )
        )
    ) / (
        1
        + 2 * gamma * pole_disp * tc
        + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
    )

    C0R = (
        (
            dSigmadr
            * (omega_tilde ^ 2 + pole_disp ^ 2)
            * r
            * (
                1
                + 2 * gamma * pole_disp * tc
                + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
            )
            * (
                -(
                    LTinv
                    * (omega_tilde ^ 2 + pole_disp ^ 2)
                    * r
                    * (1 + gamma * pole_disp * tc)
                )
                - gamma
                * LSinv
                * (omega_tilde ^ 2 + pole_disp ^ 2)
                * r
                * tc
                * (pole_disp + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc)
                - 2
                * m
                * Omega
                * omega_tilde
                * (
                    1
                    + 2 * gamma * pole_disp * tc
                    + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                )
            )
        )
        / Sigma
        - (
            -4
            * cs ^ 2
            * m ^ 2
            * Omega ^ 2
            * (-(omega_tilde ^ 2) + pole_disp ^ 2)
            * (
                1
                + 2 * gamma * pole_disp * tc
                + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
            )
            ^ 2
            * (
                gamma ^ 2 * Nr2 ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                + 2
                * D
                * gamma
                * Nr2
                * tc
                * (pole_disp + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc)
                + D ^ 2
                * (
                    1
                    + 2 * gamma * pole_disp * tc
                    + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                )
            )
            + D ^ 3
            * (
                1
                + 2 * gamma * pole_disp * tc
                + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
            )
            ^ 2
            * (
                gamma ^ 2 * omega_tilde ^ 6 * r ^ 2 * tc ^ 2
                + gamma
                * omega_tilde ^ 4
                * (
                    r ^ 2
                    + (1 + gamma) * pole_disp * r ^ 2 * tc
                    + gamma
                    * (-(cs ^ 2 * m ^ 2) + 3 * pole_disp ^ 2 * r ^ 2)
                    * tc ^ 2
                )
                + pole_disp ^ 2
                * (1 + gamma * pole_disp * tc)
                * (
                    gamma * pole_disp ^ 2 * r ^ 2 * (1 + pole_disp * tc)
                    + cs ^ 2 * m ^ 2 * (1 + gamma * pole_disp * tc)
                )
                + omega_tilde ^ 2
                * (
                    -(cs ^ 2 * m ^ 2 * (1 + 2 * gamma * pole_disp * tc))
                    + gamma
                    * pole_disp ^ 2
                    * r ^ 2
                    * (
                        2
                        + 2 * (1 + gamma) * pole_disp * tc
                        + 3 * gamma * pole_disp ^ 2 * tc ^ 2
                    )
                )
            )
            - 2
            * cs ^ 2
            * m
            * Omega
            * (
                1
                + 2 * gamma * pole_disp * tc
                + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
            )
            * (
                -2
                * gamma ^ 2
                * LSinv
                * omega_tilde
                * (omega_tilde ^ 2 + pole_disp ^ 2) ^ 2
                * r
                * tc ^ 2
                * (
                    gamma ^ 2
                    * Nr2 ^ 2
                    * (omega_tilde ^ 2 + pole_disp ^ 2)
                    * tc ^ 2
                    + 2
                    * D
                    * gamma
                    * Nr2
                    * tc
                    * (pole_disp + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc)
                    + D ^ 2
                    * (
                        1
                        + 2 * gamma * pole_disp * tc
                        + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                    )
                )
                - LTinv
                * omega_tilde
                * (omega_tilde ^ 2 + pole_disp ^ 2)
                * r
                * (1 + 2 * gamma * pole_disp * tc)
                * (
                    gamma ^ 2
                    * Nr2 ^ 2
                    * (omega_tilde ^ 2 + pole_disp ^ 2)
                    * tc ^ 2
                    + 2
                    * D
                    * gamma
                    * Nr2
                    * tc
                    * (pole_disp + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc)
                    + D ^ 2
                    * (
                        1
                        + 2 * gamma * pole_disp * tc
                        + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                    )
                )
                + r
                * (
                    gamma ^ 2
                    * Nr2 ^ 2
                    * (omega_tilde ^ 2 + pole_disp ^ 2)
                    * tc
                    * (
                        dtcdr * omega_tilde * (omega_tilde ^ 2 + pole_disp ^ 2)
                        + 2
                        * (
                            dOmegadr * m * (-(omega_tilde ^ 2) + pole_disp ^ 2)
                            + dtcdr
                            * gamma
                            * omega_tilde
                            * pole_disp
                            * (omega_tilde ^ 2 + pole_disp ^ 2)
                        )
                        * tc
                        + dOmegadr
                        * gamma
                        * m
                        * pole_disp
                        * (-5 * omega_tilde ^ 2 + 3 * pole_disp ^ 2)
                        * tc ^ 2
                        + dOmegadr
                        * gamma ^ 2
                        * m
                        * (-(omega_tilde ^ 4) + pole_disp ^ 4)
                        * tc ^ 3
                    )
                    + D
                    * (
                        1
                        + 2 * gamma * pole_disp * tc
                        + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                    )
                    * (
                        -(D * dOmegadr * gamma ^ 2 * m * omega_tilde ^ 4 * tc ^ 2)
                        + (dDdr + dNr2dr) * gamma ^ 2 * omega_tilde ^ 5 * tc ^ 2
                        + D
                        * dOmegadr
                        * m
                        * pole_disp ^ 2
                        * (1 + gamma * pole_disp * tc) ^ 2
                        - D
                        * dOmegadr
                        * m
                        * omega_tilde ^ 2
                        * (1 + 2 * gamma * pole_disp * tc)
                        + omega_tilde
                        * pole_disp ^ 2
                        * (
                            dDdr
                            + 2 * dDdr * gamma * pole_disp * tc
                            + (dDdr + dNr2dr) * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                        + omega_tilde ^ 3
                        * (
                            dDdr
                            + 2 * dDdr * gamma * pole_disp * tc
                            + 2
                            * (dDdr + dNr2dr)
                            * gamma ^ 2
                            * pole_disp ^ 2
                            * tc ^ 2
                        )
                    )
                    + gamma
                    * Nr2
                    * tc
                    * (
                        -2 * D * dOmegadr * gamma ^ 3 * m * omega_tilde ^ 6 * tc ^ 3
                        + (dDdr + dNr2dr) * gamma ^ 3 * omega_tilde ^ 7 * tc ^ 3
                        + D
                        * dOmegadr
                        * m
                        * pole_disp ^ 3
                        * (1 + gamma * pole_disp * tc) ^ 2
                        * (3 + 2 * gamma * pole_disp * tc)
                        - D
                        * dOmegadr
                        * gamma
                        * m
                        * omega_tilde ^ 4
                        * tc
                        * (
                            4
                            + 9 * gamma * pole_disp * tc
                            + 2 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                        + omega_tilde
                        * pole_disp ^ 3
                        * (1 + gamma * pole_disp * tc)
                        * (
                            2 * (dDdr + D * dtcdr * gamma * pole_disp)
                            + (3 * dDdr + dNr2dr) * gamma * pole_disp * tc
                            + (dDdr + dNr2dr) * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                        + gamma
                        * omega_tilde ^ 5
                        * (
                            2 * D * dtcdr
                            + (dDdr + dNr2dr + 2 * D * dtcdr * gamma * pole_disp) * tc
                            + 2 * (2 * dDdr + dNr2dr) * gamma * pole_disp * tc ^ 2
                            + 3
                            * (dDdr + dNr2dr)
                            * gamma ^ 2
                            * pole_disp ^ 2
                            * tc ^ 3
                        )
                        + D
                        * dOmegadr
                        * m
                        * omega_tilde ^ 2
                        * pole_disp
                        * (
                            -1
                            - 4 * gamma * pole_disp * tc
                            - 2 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                            + 2 * gamma ^ 3 * pole_disp ^ 3 * tc ^ 3
                        )
                        + omega_tilde ^ 3
                        * pole_disp
                        * (
                            2 * (dDdr + 2 * D * dtcdr * gamma * pole_disp)
                            + 2
                            * gamma
                            * pole_disp
                            * (3 * dDdr + dNr2dr + 2 * D * dtcdr * gamma * pole_disp)
                            * tc
                            + 4
                            * (2 * dDdr + dNr2dr)
                            * gamma ^ 2
                            * pole_disp ^ 2
                            * tc ^ 2
                            + 3
                            * (dDdr + dNr2dr)
                            * gamma ^ 3
                            * pole_disp ^ 3
                            * tc ^ 3
                        )
                    )
                )
            )
            + D ^ 2
            * (
                1
                + 2 * gamma * pole_disp * tc
                + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
            )
            * (
                gamma
                * Nr2
                * tc
                * (
                    3 * gamma ^ 3 * omega_tilde ^ 8 * r ^ 2 * tc ^ 3
                    + gamma
                    * omega_tilde ^ 6
                    * tc
                    * (
                        (-1 + 4 * gamma) * r ^ 2
                        + 3 * gamma * (2 + gamma) * pole_disp * r ^ 2 * tc
                        + 3
                        * gamma ^ 2
                        * (-(cs ^ 2 * m ^ 2) + 4 * pole_disp ^ 2 * r ^ 2)
                        * tc ^ 2
                    )
                    + gamma
                    * omega_tilde ^ 4
                    * (
                        3 * pole_disp * r ^ 2
                        + (
                            -3 * cs ^ 2 * m ^ 2
                            + (1 + 14 * gamma) * pole_disp ^ 2 * r ^ 2
                        )
                        * tc
                        + gamma
                        * pole_disp
                        * (
                            -7 * cs ^ 2 * m ^ 2
                            + 9 * (2 + gamma) * pole_disp ^ 2 * r ^ 2
                        )
                        * tc ^ 2
                        + 3
                        * gamma ^ 2
                        * pole_disp ^ 2
                        * (-(cs ^ 2 * m ^ 2) + 6 * pole_disp ^ 2 * r ^ 2)
                        * tc ^ 3
                    )
                    + 3
                    * pole_disp ^ 3
                    * (1 + gamma * pole_disp * tc) ^ 2
                    * (
                        gamma * pole_disp ^ 2 * r ^ 2 * (1 + pole_disp * tc)
                        + cs ^ 2 * m ^ 2 * (1 + gamma * pole_disp * tc)
                    )
                    + omega_tilde ^ 2
                    * (
                        gamma
                        * pole_disp ^ 3
                        * r ^ 2
                        * (
                            6
                            + (5 + 16 * gamma) * pole_disp * tc
                            + 9 * gamma * (2 + gamma) * pole_disp ^ 2 * tc ^ 2
                            + 12 * gamma ^ 2 * pole_disp ^ 3 * tc ^ 3
                        )
                        + cs ^ 2
                        * m ^ 2
                        * pole_disp
                        * (
                            -1
                            - 2 * gamma * pole_disp * tc
                            + 2 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                            + 3 * gamma ^ 3 * pole_disp ^ 3 * tc ^ 3
                        )
                    )
                )
                + cs ^ 2
                * (omega_tilde ^ 2 + pole_disp ^ 2)
                * (
                    -(
                        r
                        * (
                            1
                            + 2 * gamma * pole_disp * tc
                            + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                        )
                        * (
                            -2 * dOmegadr * m * omega_tilde
                            + gamma
                            * pole_disp
                            * (
                                -4 * dOmegadr * m * omega_tilde
                                + (omega_tilde ^ 2 + pole_disp ^ 2)
                                * (-LSinv - dLSinvdr * r)
                            )
                            * tc
                            + gamma ^ 2
                            * (omega_tilde ^ 2 + pole_disp ^ 2)
                            * (
                                -2 * dOmegadr * m * omega_tilde
                                + (omega_tilde ^ 2 + pole_disp ^ 2)
                                * (-LSinv - dLSinvdr * r)
                            )
                            * tc ^ 2
                            + (omega_tilde ^ 2 + pole_disp ^ 2)
                            * (-LTinv - dLTinvdr * r)
                            * (1 + gamma * pole_disp * tc)
                        )
                    )
                    + gamma
                    * LSinv
                    * (omega_tilde ^ 2 + pole_disp ^ 2)
                    * r ^ 2
                    * (
                        -2
                        * dOmegadr
                        * gamma
                        * m
                        * omega_tilde
                        * tc ^ 2
                        * (1 + gamma * pole_disp * tc)
                        + dtcdr * pole_disp * (1 + gamma * pole_disp * tc) ^ 2
                        + dtcdr
                        * gamma
                        * omega_tilde ^ 2
                        * tc
                        * (2 + gamma * pole_disp * tc)
                    )
                    + gamma ^ 2
                    * LSinv ^ 2
                    * (omega_tilde ^ 2 + pole_disp ^ 2)
                    * r ^ 2
                    * tc ^ 2
                    * (
                        gamma ^ 2 * omega_tilde ^ 4 * tc ^ 2
                        + pole_disp ^ 2 * (1 + gamma * pole_disp * tc) ^ 2
                        + omega_tilde ^ 2
                        * (
                            -1
                            + 2 * gamma * pole_disp * tc
                            + 2 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                    )
                    - gamma
                    * LTinv
                    * (omega_tilde ^ 2 + pole_disp ^ 2)
                    * r
                    * (
                        -(
                            LSinv
                            * r
                            * tc
                            * (
                                pole_disp
                                + 2 * gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc
                                + gamma ^ 2
                                * pole_disp
                                * (omega_tilde ^ 2 + pole_disp ^ 2)
                                * tc ^ 2
                            )
                        )
                        + r
                        * (
                            -2
                            * dOmegadr
                            * gamma
                            * m
                            * omega_tilde
                            * tc ^ 2
                            * (1 + gamma * pole_disp * tc)
                            + dtcdr * pole_disp * (1 + gamma * pole_disp * tc) ^ 2
                            + dtcdr
                            * gamma
                            * omega_tilde ^ 2
                            * tc
                            * (2 + gamma * pole_disp * tc)
                        )
                    )
                )
            )
            + gamma
            * Nr2
            * (omega_tilde ^ 2 + pole_disp ^ 2) ^ 2
            * tc
            * (
                cs ^ 2
                * r
                * (
                    1
                    + 2 * gamma * pole_disp * tc
                    + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                )
                * (
                    -(
                        gamma
                        * LSinv
                        * (omega_tilde ^ 2 + pole_disp ^ 2)
                        * r
                        * tc
                        * (
                            dDdr
                            + (2 * dDdr + dNr2dr) * gamma * pole_disp * tc
                            + (dDdr + dNr2dr)
                            * gamma ^ 2
                            * (omega_tilde ^ 2 + pole_disp ^ 2)
                            * tc ^ 2
                        )
                    )
                    - LTinv
                    * r
                    * (
                        dDdr * pole_disp
                        + gamma
                        * (
                            2 * dDdr * pole_disp ^ 2
                            + dNr2dr * (omega_tilde ^ 2 + pole_disp ^ 2)
                        )
                        * tc
                        + (dDdr + dNr2dr)
                        * gamma ^ 2
                        * pole_disp
                        * (omega_tilde ^ 2 + pole_disp ^ 2)
                        * tc ^ 2
                    )
                )
                + gamma ^ 2
                * Nr2 ^ 2
                * tc ^ 2
                * (
                    gamma ^ 3 * omega_tilde ^ 6 * r ^ 2 * tc ^ 3
                    + gamma
                    * omega_tilde ^ 4
                    * tc
                    * (
                        (-1 + 2 * gamma) * r ^ 2
                        + gamma * (2 + gamma) * pole_disp * r ^ 2 * tc
                        + gamma ^ 2
                        * (-(cs ^ 2 * m ^ 2) + 3 * pole_disp ^ 2 * r ^ 2)
                        * tc ^ 2
                    )
                    + gamma
                    * omega_tilde ^ 2
                    * (
                        pole_disp * r ^ 2
                        + (-(cs ^ 2 * m ^ 2) + 4 * gamma * pole_disp ^ 2 * r ^ 2)
                        * tc
                        + gamma
                        * pole_disp
                        * (
                            -(cs ^ 2 * m ^ 2)
                            + 2 * (2 + gamma) * pole_disp ^ 2 * r ^ 2
                        )
                        * tc ^ 2
                        + 3 * gamma ^ 2 * pole_disp ^ 4 * r ^ 2 * tc ^ 3
                    )
                    + pole_disp
                    * (1 + gamma * pole_disp * tc) ^ 2
                    * (
                        gamma * pole_disp ^ 2 * r ^ 2 * (1 + pole_disp * tc)
                        + cs ^ 2 * m ^ 2 * (1 + gamma * pole_disp * tc)
                    )
                )
                + cs ^ 2
                * gamma
                * Nr2
                * (
                    -(
                        r
                        * tc
                        * (
                            1
                            + 2 * gamma * pole_disp * tc
                            + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                        )
                        * (
                            -2 * dOmegadr * m * omega_tilde
                            + gamma
                            * pole_disp
                            * (
                                -4 * dOmegadr * m * omega_tilde
                                + (omega_tilde ^ 2 + pole_disp ^ 2)
                                * (-LSinv - dLSinvdr * r)
                            )
                            * tc
                            + gamma ^ 2
                            * (omega_tilde ^ 2 + pole_disp ^ 2)
                            * (
                                -2 * dOmegadr * m * omega_tilde
                                + (omega_tilde ^ 2 + pole_disp ^ 2)
                                * (-LSinv - dLSinvdr * r)
                            )
                            * tc ^ 2
                            + (omega_tilde ^ 2 + pole_disp ^ 2)
                            * (-LTinv - dLTinvdr * r)
                            * (1 + gamma * pole_disp * tc)
                        )
                    )
                    + gamma ^ 2
                    * LSinv ^ 2
                    * (omega_tilde ^ 2 + pole_disp ^ 2)
                    * r ^ 2
                    * tc ^ 3
                    * (
                        gamma ^ 2 * omega_tilde ^ 4 * tc ^ 2
                        + pole_disp ^ 2 * (1 + gamma * pole_disp * tc) ^ 2
                        + omega_tilde ^ 2
                        * (
                            -1
                            + 2 * gamma * pole_disp * tc
                            + 2 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                    )
                    - LTinv
                    * r
                    * (
                        -(
                            gamma
                            * LSinv
                            * (omega_tilde ^ 2 + pole_disp ^ 2)
                            * r
                            * tc ^ 2
                            * (
                                pole_disp
                                + 2 * gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc
                                + gamma ^ 2
                                * pole_disp
                                * (omega_tilde ^ 2 + pole_disp ^ 2)
                                * tc ^ 2
                            )
                        )
                        + r
                        * (
                            1
                            + 2 * gamma * pole_disp * tc
                            + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                        )
                        * (
                            dtcdr * omega_tilde ^ 2 * (1 + gamma * pole_disp * tc)
                            + dtcdr * pole_disp ^ 2 * (1 + gamma * pole_disp * tc)
                            - dOmegadr
                            * m
                            * omega_tilde
                            * tc
                            * (1 + 2 * gamma * pole_disp * tc)
                        )
                    )
                )
            )
            + D
            * (omega_tilde ^ 2 + pole_disp ^ 2)
            * (
                cs ^ 2
                * (omega_tilde ^ 2 + pole_disp ^ 2)
                * r
                * (
                    1
                    + 2 * gamma * pole_disp * tc
                    + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                )
                * (
                    -(
                        LTinv
                        * r
                        * (
                            dDdr
                            + (3 * dDdr + dNr2dr) * gamma * pole_disp * tc
                            + gamma ^ 2
                            * (
                                2 * dNr2dr * (omega_tilde ^ 2 + pole_disp ^ 2)
                                + dDdr * (omega_tilde ^ 2 + 3 * pole_disp ^ 2)
                            )
                            * tc ^ 2
                            + (dDdr + dNr2dr)
                            * gamma ^ 3
                            * pole_disp
                            * (omega_tilde ^ 2 + pole_disp ^ 2)
                            * tc ^ 3
                        )
                    )
                    - gamma
                    * LSinv
                    * r
                    * tc
                    * (
                        dDdr * pole_disp
                        + gamma
                        * (
                            dNr2dr * (-(omega_tilde ^ 2) + pole_disp ^ 2)
                            + dDdr * (omega_tilde ^ 2 + 3 * pole_disp ^ 2)
                        )
                        * tc
                        + (3 * dDdr + 2 * dNr2dr)
                        * gamma ^ 2
                        * pole_disp
                        * (omega_tilde ^ 2 + pole_disp ^ 2)
                        * tc ^ 2
                        + (dDdr + dNr2dr)
                        * gamma ^ 3
                        * (omega_tilde ^ 2 + pole_disp ^ 2) ^ 2
                        * tc ^ 3
                    )
                )
                + gamma ^ 2
                * Nr2 ^ 2
                * tc ^ 2
                * (
                    3 * gamma ^ 4 * omega_tilde ^ 8 * r ^ 2 * tc ^ 4
                    + gamma ^ 2
                    * omega_tilde ^ 6
                    * tc ^ 2
                    * (
                        (-1 + 5 * gamma) * r ^ 2
                        + 3 * gamma * (3 + gamma) * pole_disp * r ^ 2 * tc
                        + 3
                        * gamma ^ 2
                        * (-(cs ^ 2 * m ^ 2) + 4 * pole_disp ^ 2 * r ^ 2)
                        * tc ^ 2
                    )
                    + gamma
                    * omega_tilde ^ 4
                    * (
                        r ^ 2
                        + (-1 + 9 * gamma) * pole_disp * r ^ 2 * tc
                        + gamma
                        * (
                            -4 * cs ^ 2 * m ^ 2
                            + (7 + 19 * gamma) * pole_disp ^ 2 * r ^ 2
                        )
                        * tc ^ 2
                        + gamma ^ 2
                        * pole_disp
                        * (
                            -8 * cs ^ 2 * m ^ 2
                            + 9 * (3 + gamma) * pole_disp ^ 2 * r ^ 2
                        )
                        * tc ^ 3
                        + 3
                        * gamma ^ 3
                        * pole_disp ^ 2
                        * (-(cs ^ 2 * m ^ 2) + 6 * pole_disp ^ 2 * r ^ 2)
                        * tc ^ 4
                    )
                    + 3
                    * pole_disp ^ 2
                    * (1 + gamma * pole_disp * tc) ^ 3
                    * (
                        gamma * pole_disp ^ 2 * r ^ 2 * (1 + pole_disp * tc)
                        + cs ^ 2 * m ^ 2 * (1 + gamma * pole_disp * tc)
                    )
                    + omega_tilde ^ 2
                    * (1 + gamma * pole_disp * tc)
                    * (
                        gamma
                        * pole_disp ^ 2
                        * r ^ 2
                        * (
                            4
                            + 2 * (pole_disp + 7 * gamma * pole_disp) * tc
                            + 3 * gamma * (5 + 3 * gamma) * pole_disp ^ 2 * tc ^ 2
                            + 12 * gamma ^ 2 * pole_disp ^ 3 * tc ^ 3
                        )
                        + cs ^ 2
                        * m ^ 2
                        * (
                            -1
                            - 3 * gamma * pole_disp * tc
                            + gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                            + 3 * gamma ^ 3 * pole_disp ^ 3 * tc ^ 3
                        )
                    )
                )
                + cs ^ 2
                * gamma
                * Nr2
                * (
                    -(
                        LTinv
                        * (omega_tilde ^ 2 + pole_disp ^ 2)
                        * r
                        * (
                            -2
                            * gamma
                            * LSinv
                            * r
                            * tc ^ 2
                            * (
                                pole_disp ^ 2
                                + 3
                                * gamma
                                * pole_disp
                                * (omega_tilde ^ 2 + pole_disp ^ 2)
                                * tc
                                + gamma ^ 2
                                * (
                                    2 * omega_tilde ^ 4
                                    + 5 * omega_tilde ^ 2 * pole_disp ^ 2
                                    + 3 * pole_disp ^ 4
                                )
                                * tc ^ 2
                                + gamma ^ 3
                                * pole_disp
                                * (omega_tilde ^ 2 + pole_disp ^ 2) ^ 2
                                * tc ^ 3
                            )
                            + r
                            * (
                                1
                                + 2 * gamma * pole_disp * tc
                                + gamma ^ 2
                                * (omega_tilde ^ 2 + pole_disp ^ 2)
                                * tc ^ 2
                            )
                            * (
                                dtcdr
                                * gamma
                                * omega_tilde ^ 2
                                * tc
                                * (3 + 2 * gamma * pole_disp * tc)
                                - dOmegadr
                                * gamma
                                * m
                                * omega_tilde
                                * tc ^ 2
                                * (3 + 4 * gamma * pole_disp * tc)
                                + dtcdr
                                * pole_disp
                                * (
                                    1
                                    + 3 * gamma * pole_disp * tc
                                    + 2 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                                )
                            )
                        )
                    )
                    + tc
                    * (
                        -2
                        * r
                        * (
                            pole_disp
                            + gamma * (omega_tilde ^ 2 + 3 * pole_disp ^ 2) * tc
                            + 3
                            * gamma ^ 2
                            * pole_disp
                            * (omega_tilde ^ 2 + pole_disp ^ 2)
                            * tc ^ 2
                            + gamma ^ 3
                            * (omega_tilde ^ 2 + pole_disp ^ 2) ^ 2
                            * tc ^ 3
                        )
                        * (
                            -2 * dOmegadr * m * omega_tilde
                            + gamma
                            * pole_disp
                            * (
                                -4 * dOmegadr * m * omega_tilde
                                + (omega_tilde ^ 2 + pole_disp ^ 2)
                                * (-LSinv - dLSinvdr * r)
                            )
                            * tc
                            + gamma ^ 2
                            * (omega_tilde ^ 2 + pole_disp ^ 2)
                            * (
                                -2 * dOmegadr * m * omega_tilde
                                + (omega_tilde ^ 2 + pole_disp ^ 2)
                                * (-LSinv - dLSinvdr * r)
                            )
                            * tc ^ 2
                            + (omega_tilde ^ 2 + pole_disp ^ 2)
                            * (-LTinv - dLTinvdr * r)
                            * (1 + gamma * pole_disp * tc)
                        )
                        + gamma
                        * LSinv
                        * (omega_tilde ^ 2 + pole_disp ^ 2)
                        * r ^ 2
                        * (
                            1
                            + 2 * gamma * pole_disp * tc
                            + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                        )
                        * (
                            dtcdr * omega_tilde ^ 2 * (1 + gamma * pole_disp * tc)
                            + dtcdr * pole_disp ^ 2 * (1 + gamma * pole_disp * tc)
                            - dOmegadr
                            * m
                            * omega_tilde
                            * tc
                            * (1 + 2 * gamma * pole_disp * tc)
                        )
                        + 2
                        * gamma ^ 2
                        * LSinv ^ 2
                        * (omega_tilde ^ 2 + pole_disp ^ 2)
                        * r ^ 2
                        * tc ^ 2
                        * (
                            gamma ^ 3 * omega_tilde ^ 6 * tc ^ 3
                            + pole_disp ^ 3 * (1 + gamma * pole_disp * tc) ^ 3
                            + gamma
                            * omega_tilde ^ 4
                            * tc
                            * (
                                -1
                                + 3 * gamma * pole_disp * tc
                                + 3 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                            )
                            + omega_tilde ^ 2
                            * pole_disp
                            * (
                                -1
                                + 2 * gamma * pole_disp * tc
                                + 6 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                                + 3 * gamma ^ 3 * pole_disp ^ 3 * tc ^ 3
                            )
                        )
                    )
                )
            )
        )
        / (
            cs ^ 2
            * (
                gamma ^ 2 * Nr2 ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                + 2
                * D
                * gamma
                * Nr2
                * tc
                * (pole_disp + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc)
                + D ^ 2
                * (
                    1
                    + 2 * gamma * pole_disp * tc
                    + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                )
            )
        )
    ) / (
        (omega_tilde ^ 2 + pole_disp ^ 2) ^ 2
        * r ^ 2
        * (
            1
            + 2 * gamma * pole_disp * tc
            + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
        )
        ^ 2
    )

    C0I = (
        cs ^ 2
        * dSigmadr
        * (omega_tilde ^ 2 + pole_disp ^ 2)
        * r
        * (
            1
            + 2 * gamma * pole_disp * tc
            + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
        )
        * (
            gamma ^ 2 * Nr2 ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
            + 2
            * D
            * gamma
            * Nr2
            * tc
            * (pole_disp + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc)
            + D ^ 2
            * (
                1
                + 2 * gamma * pole_disp * tc
                + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
            )
        )
        * (
            gamma * LSinv * omega_tilde * (omega_tilde ^ 2 + pole_disp ^ 2) * r * tc
            + 2
            * m
            * Omega
            * pole_disp
            * (
                1
                + 2 * gamma * pole_disp * tc
                + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
            )
        )
        + Sigma
        * (
            cs ^ 2
            * dOmegadr
            * gamma ^ 5
            * m
            * (D + Nr2)
            * omega_tilde ^ 8
            * r
            * tc ^ 5
            * (D * LSinv * r + 2 * gamma * (D + Nr2) * pole_disp * tc)
            + cs ^ 2
            * dOmegadr
            * m
            * pole_disp ^ 3
            * r
            * (1 + gamma * pole_disp * tc) ^ 3
            * (gamma * Nr2 * pole_disp * tc + D * (1 + gamma * pole_disp * tc))
            * (
                2 * gamma * Nr2 * pole_disp * tc * (1 + gamma * pole_disp * tc)
                + D
                * (
                    2
                    + gamma * pole_disp * (4 - LSinv * r) * tc
                    + 2 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                )
            )
            + 8
            * cs ^ 2
            * m ^ 2
            * Omega ^ 2
            * omega_tilde
            * pole_disp
            * (
                1
                + 2 * gamma * pole_disp * tc
                + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
            )
            ^ 2
            * (
                gamma ^ 2 * Nr2 ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                + 2
                * D
                * gamma
                * Nr2
                * tc
                * (pole_disp + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc)
                + D ^ 2
                * (
                    1
                    + 2 * gamma * pole_disp * tc
                    + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                )
            )
            + cs ^ 2
            * dOmegadr
            * m
            * omega_tilde ^ 2
            * pole_disp
            * r
            * (1 + gamma * pole_disp * tc)
            * (
                2
                * D ^ 2
                * (1 + gamma * pole_disp * tc) ^ 3
                * (
                    1
                    + gamma * pole_disp * (2 - LSinv * r) * tc
                    + 4 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                )
                + 4
                * gamma ^ 2
                * Nr2 ^ 2
                * pole_disp ^ 2
                * tc ^ 2
                * (
                    1
                    + 3 * gamma * pole_disp * tc
                    + 4 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                    + 2 * gamma ^ 3 * pole_disp ^ 3 * tc ^ 3
                )
                + D
                * gamma
                * Nr2
                * pole_disp
                * tc
                * (
                    4
                    + 2 * gamma * pole_disp * (10 - LSinv * r) * tc
                    + gamma ^ 2 * pole_disp ^ 2 * (44 - 3 * LSinv * r) * tc ^ 2
                    + 2 * gamma ^ 3 * pole_disp ^ 3 * (22 - LSinv * r) * tc ^ 3
                    + 16 * gamma ^ 4 * pole_disp ^ 4 * tc ^ 4
                )
            )
            + gamma ^ 5
            * omega_tilde ^ 9
            * tc ^ 4
            * (
                -((-2 + gamma) * Nr2 ^ 3 * r ^ 2 * tc)
                + Nr2 ^ 2
                * (
                    2 * cs ^ 2 * LSinv ^ 2 * r ^ 2
                    + r
                    * (
                        D * (5 - 3 * gamma) * r
                        + cs ^ 2 * (LSinv - LTinv + dLSinvdr * r - dLTinvdr * r)
                    )
                )
                * tc
                + Nr2
                * (
                    4 * cs ^ 2 * D * LSinv ^ 2 * r ^ 2 * tc
                    + D
                    * r
                    * (
                        D * (4 - 3 * gamma) * r
                        + 2 * cs ^ 2 * (LSinv - LTinv + dLSinvdr * r - dLTinvdr * r)
                    )
                    * tc
                    - cs ^ 2 * LSinv * r ^ 2 * (D * dtcdr + dNr2dr * tc)
                )
                + D
                * (
                    2 * cs ^ 2 * D * LSinv ^ 2 * r ^ 2 * tc
                    + D
                    * r
                    * (
                        cs ^ 2 * (LSinv - LTinv + dLSinvdr * r - dLTinvdr * r)
                        + D * (r - gamma * r)
                    )
                    * tc
                    - cs ^ 2 * LSinv * r ^ 2 * (D * dtcdr + (dDdr + 2 * dNr2dr) * tc)
                )
            )
            + cs ^ 2
            * dOmegadr
            * gamma ^ 3
            * m
            * omega_tilde ^ 6
            * r
            * tc ^ 3
            * (
                4
                * gamma
                * Nr2 ^ 2
                * pole_disp
                * tc
                * (
                    1
                    + 2 * gamma * pole_disp * tc
                    + 2 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                )
                + 2
                * D ^ 2
                * gamma
                * pole_disp
                * tc
                * (
                    3
                    - gamma * pole_disp * (-6 - LSinv * r) * tc
                    + 4 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                )
                - D
                * Nr2
                * (
                    -(
                        LSinv
                        * r
                        * (
                            1
                            + gamma * pole_disp * tc
                            + 2 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                    )
                    - 4
                    * gamma
                    * pole_disp
                    * tc
                    * (
                        2
                        + 5 * gamma * pole_disp * tc
                        + 4 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                    )
                )
            )
            + cs ^ 2
            * dOmegadr
            * gamma
            * m
            * omega_tilde ^ 4
            * r
            * tc
            * (
                -(
                    D
                    * LSinv
                    * r
                    * (
                        gamma
                        * Nr2
                        * pole_disp
                        * tc
                        * (
                            1
                            + gamma * pole_disp * tc
                            + gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                        + D
                        * (
                            1
                            + 4 * gamma * pole_disp * tc
                            + 6 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                            + 4 * gamma ^ 3 * pole_disp ^ 3 * tc ^ 3
                        )
                    )
                )
                + 2
                * gamma
                * pole_disp
                * tc
                * (
                    3
                    * D ^ 2
                    * (1 + gamma * pole_disp * tc) ^ 2
                    * (
                        1
                        + 2 * gamma * pole_disp * tc
                        + 2 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                    )
                    + Nr2 ^ 2
                    * (
                        1
                        + 4 * gamma * pole_disp * tc
                        + 10 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        + 12 * gamma ^ 3 * pole_disp ^ 3 * tc ^ 3
                        + 6 * gamma ^ 4 * pole_disp ^ 4 * tc ^ 4
                    )
                    + 2
                    * D
                    * Nr2
                    * (
                        1
                        + 6 * gamma * pole_disp * tc
                        + 14 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        + 15 * gamma ^ 3 * pole_disp ^ 3 * tc ^ 3
                        + 6 * gamma ^ 4 * pole_disp ^ 4 * tc ^ 4
                    )
                )
            )
            - 2
            * cs ^ 2
            * m
            * Omega
            * (
                1
                + 2 * gamma * pole_disp * tc
                + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
            )
            * (
                -2
                * gamma
                * LSinv
                * (omega_tilde ^ 2 + pole_disp ^ 2) ^ 2
                * r
                * tc
                * (1 + gamma * pole_disp * tc)
                * (
                    gamma ^ 2
                    * Nr2 ^ 2
                    * (omega_tilde ^ 2 + pole_disp ^ 2)
                    * tc ^ 2
                    + 2
                    * D
                    * gamma
                    * Nr2
                    * tc
                    * (pole_disp + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc)
                    + D ^ 2
                    * (
                        1
                        + 2 * gamma * pole_disp * tc
                        + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                    )
                )
                + r
                * (
                    -2
                    * D ^ 2
                    * dOmegadr
                    * m
                    * omega_tilde
                    * pole_disp
                    * (
                        1
                        + 2 * gamma * pole_disp * tc
                        + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                    )
                    ^ 2
                    + gamma
                    * Nr2
                    * (omega_tilde ^ 2 + pole_disp ^ 2)
                    * tc
                    * (
                        -(dDdr * omega_tilde ^ 2)
                        + dtcdr * gamma * Nr2 * omega_tilde ^ 2 * pole_disp
                        + dDdr * pole_disp ^ 2
                        + dtcdr * gamma * Nr2 * pole_disp ^ 3
                        + gamma
                        * (
                            dNr2dr * pole_disp * (omega_tilde ^ 2 + pole_disp ^ 2)
                            + dDdr
                            * pole_disp
                            * (-(omega_tilde ^ 2) + 3 * pole_disp ^ 2)
                            + Nr2
                            * (
                                -(dtcdr * gamma * omega_tilde ^ 4)
                                - 4 * dOmegadr * m * omega_tilde * pole_disp
                                + dtcdr * gamma * pole_disp ^ 4
                            )
                        )
                        * tc
                        + gamma ^ 2
                        * (
                            dOmegadr
                            * m
                            * Nr2
                            * omega_tilde
                            * (omega_tilde ^ 2 - 7 * pole_disp ^ 2)
                            + 2
                            * dNr2dr
                            * pole_disp ^ 2
                            * (omega_tilde ^ 2 + pole_disp ^ 2)
                            + dDdr
                            * (
                                -(omega_tilde ^ 4)
                                + 2 * omega_tilde ^ 2 * pole_disp ^ 2
                                + 3 * pole_disp ^ 4
                            )
                        )
                        * tc ^ 2
                        + gamma ^ 3
                        * pole_disp
                        * (omega_tilde ^ 2 + pole_disp ^ 2)
                        * (
                            -2 * dOmegadr * m * Nr2 * omega_tilde
                            + dDdr * (omega_tilde ^ 2 + pole_disp ^ 2)
                            + dNr2dr * (omega_tilde ^ 2 + pole_disp ^ 2)
                        )
                        * tc ^ 3
                    )
                    + D
                    * (
                        (omega_tilde ^ 2 + pole_disp ^ 2)
                        * (
                            dDdr * pole_disp
                            + dtcdr * gamma * Nr2 * (omega_tilde ^ 2 + pole_disp ^ 2)
                        )
                        + gamma
                        * (
                            4
                            * dDdr
                            * pole_disp ^ 2
                            * (omega_tilde ^ 2 + pole_disp ^ 2)
                            + dNr2dr * (omega_tilde ^ 2 + pole_disp ^ 2) ^ 2
                            + Nr2
                            * (
                                -(dOmegadr * m * omega_tilde ^ 3)
                                + 2 * dtcdr * gamma * omega_tilde ^ 4 * pole_disp
                                - 5 * dOmegadr * m * omega_tilde * pole_disp ^ 2
                                + 4 * dtcdr * gamma * omega_tilde ^ 2 * pole_disp ^ 3
                                + 2 * dtcdr * gamma * pole_disp ^ 5
                            )
                        )
                        * tc
                        + gamma ^ 2
                        * (
                            3
                            * dNr2dr
                            * pole_disp
                            * (omega_tilde ^ 2 + pole_disp ^ 2) ^ 2
                            + 2
                            * dDdr
                            * pole_disp
                            * (
                                omega_tilde ^ 4
                                + 4 * omega_tilde ^ 2 * pole_disp ^ 2
                                + 3 * pole_disp ^ 4
                            )
                            - Nr2
                            * (
                                dtcdr * gamma * omega_tilde ^ 6
                                + 8 * dOmegadr * m * omega_tilde ^ 3 * pole_disp
                                + dtcdr * gamma * omega_tilde ^ 4 * pole_disp ^ 2
                                + 16 * dOmegadr * m * omega_tilde * pole_disp ^ 3
                                - dtcdr * gamma * omega_tilde ^ 2 * pole_disp ^ 4
                                - dtcdr * gamma * pole_disp ^ 6
                            )
                        )
                        * tc ^ 2
                        + gamma ^ 3
                        * (omega_tilde ^ 2 + pole_disp ^ 2)
                        * (
                            dOmegadr
                            * m
                            * Nr2
                            * omega_tilde
                            * (omega_tilde ^ 2 - 15 * pole_disp ^ 2)
                            + 4
                            * dDdr
                            * pole_disp ^ 2
                            * (omega_tilde ^ 2 + pole_disp ^ 2)
                            + dNr2dr
                            * (
                                omega_tilde ^ 4
                                + 4 * omega_tilde ^ 2 * pole_disp ^ 2
                                + 3 * pole_disp ^ 4
                            )
                        )
                        * tc ^ 3
                        + gamma ^ 4
                        * pole_disp
                        * (omega_tilde ^ 2 + pole_disp ^ 2) ^ 2
                        * (
                            -4 * dOmegadr * m * Nr2 * omega_tilde
                            + dDdr * (omega_tilde ^ 2 + pole_disp ^ 2)
                            + dNr2dr * (omega_tilde ^ 2 + pole_disp ^ 2)
                        )
                        * tc ^ 4
                    )
                )
            )
            + gamma ^ 3
            * omega_tilde ^ 7
            * tc ^ 2
            * (
                gamma
                * Nr2 ^ 3
                * tc
                * (
                    r ^ 2
                    + 2 * pole_disp * r ^ 2 * tc
                    - gamma
                    * (cs ^ 2 * m ^ 2 + 4 * (-2 + gamma) * pole_disp ^ 2 * r ^ 2)
                    * tc ^ 2
                    - 2 * cs ^ 2 * gamma ^ 2 * m ^ 2 * pole_disp * tc ^ 3
                )
                + Nr2 ^ 2
                * tc
                * (
                    D
                    * (
                        (1 + gamma) * r ^ 2
                        + 2 * (5 - 2 * gamma) * gamma * pole_disp * r ^ 2 * tc
                        - 2
                        * gamma ^ 2
                        * (
                            cs ^ 2 * m ^ 2
                            + 2 * (-5 + 3 * gamma) * pole_disp ^ 2 * r ^ 2
                        )
                        * tc ^ 2
                        - 6 * cs ^ 2 * gamma ^ 3 * m ^ 2 * pole_disp * tc ^ 3
                    )
                    + cs ^ 2
                    * (
                        2
                        * gamma
                        * LSinv ^ 2
                        * pole_disp
                        * r ^ 2
                        * tc
                        * (1 + 4 * gamma * pole_disp * tc)
                        + r
                        * (LSinv - LTinv + dLSinvdr * r - dLTinvdr * r)
                        * (
                            1
                            + 2 * gamma * pole_disp * tc
                            + 4 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                    )
                )
                + Nr2
                * (
                    -(
                        cs ^ 2
                        * dNr2dr
                        * LSinv
                        * r ^ 2
                        * tc
                        * (
                            1
                            + 2 * gamma * pole_disp * tc
                            + 4 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                    )
                    - D ^ 2
                    * tc
                    * (
                        2 * (-2 + gamma) * r ^ 2
                        + 4 * gamma * (-3 + 2 * gamma) * pole_disp * r ^ 2 * tc
                        + gamma ^ 2
                        * (
                            cs ^ 2 * m ^ 2
                            + 4 * (-4 + 3 * gamma) * pole_disp ^ 2 * r ^ 2
                        )
                        * tc ^ 2
                        + 6 * cs ^ 2 * gamma ^ 3 * m ^ 2 * pole_disp * tc ^ 3
                    )
                    + cs ^ 2
                    * D
                    * (
                        8
                        * gamma
                        * LSinv ^ 2
                        * pole_disp
                        * r ^ 2
                        * tc ^ 2
                        * (1 + 2 * gamma * pole_disp * tc)
                        - dtcdr
                        * LSinv
                        * r ^ 2
                        * (
                            1
                            + 2 * gamma * pole_disp * tc
                            + 4 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                        + 2
                        * r
                        * (LSinv - LTinv + dLSinvdr * r - dLTinvdr * r)
                        * tc
                        * (
                            1
                            + 3 * gamma * pole_disp * tc
                            + 4 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                    )
                )
                + 2
                * D
                * tc
                * (
                    cs ^ 2
                    * D
                    * LSinv ^ 2
                    * r ^ 2
                    * (
                        1
                        + 3 * gamma * pole_disp * tc
                        + 4 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                    )
                    - cs ^ 2
                    * LSinv
                    * r ^ 2
                    * (
                        D * dtcdr * gamma * pole_disp * (1 + 2 * gamma * pole_disp * tc)
                        + dDdr
                        * (
                            1
                            + 2 * gamma * pole_disp * tc
                            + 2 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                        + dNr2dr
                        * (
                            1
                            + 3 * gamma * pole_disp * tc
                            + 4 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                    )
                    + D
                    * (
                        cs ^ 2
                        * r
                        * (LSinv - LTinv + dLSinvdr * r - dLTinvdr * r)
                        * (
                            1
                            + 2 * gamma * pole_disp * tc
                            + 2 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                        - D
                        * (
                            (-1 + gamma) * r ^ 2
                            + 2 * (-1 + gamma) * gamma * pole_disp * r ^ 2 * tc
                            + 2
                            * (-1 + gamma)
                            * gamma ^ 2
                            * pole_disp ^ 2
                            * r ^ 2
                            * tc ^ 2
                            + cs ^ 2 * gamma ^ 3 * m ^ 2 * pole_disp * tc ^ 3
                        )
                    )
                )
            )
            + gamma
            * omega_tilde ^ 5
            * (
                -(
                    gamma ^ 2
                    * Nr2 ^ 3
                    * tc ^ 3
                    * (
                        3
                        * gamma
                        * pole_disp ^ 2
                        * r ^ 2
                        * (
                            -1
                            - 2 * pole_disp * tc
                            + 2 * (-2 + gamma) * gamma * pole_disp ^ 2 * tc ^ 2
                        )
                        + cs ^ 2
                        * m ^ 2
                        * (
                            1
                            + 4 * gamma * pole_disp * tc
                            + 7 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                            + 6 * gamma ^ 3 * pole_disp ^ 3 * tc ^ 3
                        )
                    )
                )
                + gamma ^ 2
                * Nr2 ^ 2
                * tc ^ 2
                * (
                    D
                    * (
                        2 * pole_disp * r ^ 2
                        + (
                            -2 * cs ^ 2 * m ^ 2
                            + (7 + 3 * gamma) * pole_disp ^ 2 * r ^ 2
                        )
                        * tc
                        - 2
                        * gamma
                        * pole_disp
                        * (
                            7 * cs ^ 2 * m ^ 2
                            + 3 * (-5 + 2 * gamma) * pole_disp ^ 2 * r ^ 2
                        )
                        * tc ^ 2
                        - 2
                        * gamma ^ 2
                        * pole_disp ^ 2
                        * (
                            13 * cs ^ 2 * m ^ 2
                            + 3 * (-5 + 3 * gamma) * pole_disp ^ 2 * r ^ 2
                        )
                        * tc ^ 3
                        - 18 * cs ^ 2 * gamma ^ 3 * m ^ 2 * pole_disp ^ 3 * tc ^ 4
                    )
                    + 3
                    * cs ^ 2
                    * pole_disp ^ 2
                    * tc
                    * (
                        2
                        * gamma
                        * LSinv ^ 2
                        * pole_disp
                        * r ^ 2
                        * tc
                        * (1 + 2 * gamma * pole_disp * tc)
                        + r
                        * (LSinv - LTinv + dLSinvdr * r - dLTinvdr * r)
                        * (
                            1
                            + 2 * gamma * pole_disp * tc
                            + 2 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                    )
                )
                + gamma
                * Nr2
                * tc
                * (
                    -3
                    * cs ^ 2
                    * dNr2dr
                    * gamma
                    * LSinv
                    * pole_disp ^ 2
                    * r ^ 2
                    * tc ^ 2
                    * (
                        1
                        + 2 * gamma * pole_disp * tc
                        + 2 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                    )
                    + D ^ 2
                    * (
                        r ^ 2
                        + 4 * pole_disp * r ^ 2 * tc
                        - 2
                        * gamma
                        * (
                            cs ^ 2 * m ^ 2
                            + 5 * (-2 + gamma) * pole_disp ^ 2 * r ^ 2
                        )
                        * tc ^ 2
                        - 4
                        * gamma ^ 2
                        * pole_disp
                        * (
                            4 * cs ^ 2 * m ^ 2
                            + 3 * (-3 + 2 * gamma) * pole_disp ^ 2 * r ^ 2
                        )
                        * tc ^ 3
                        + gamma ^ 3
                        * pole_disp ^ 2
                        * (
                            -31 * cs ^ 2 * m ^ 2
                            + 6 * (4 - 3 * gamma) * pole_disp ^ 2 * r ^ 2
                        )
                        * tc ^ 4
                        - 18 * cs ^ 2 * gamma ^ 4 * m ^ 2 * pole_disp ^ 3 * tc ^ 5
                    )
                    + cs ^ 2
                    * D
                    * pole_disp
                    * tc
                    * (
                        -3
                        * dtcdr
                        * gamma
                        * LSinv
                        * pole_disp
                        * r ^ 2
                        * (
                            1
                            + 2 * gamma * pole_disp * tc
                            + 2 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                        + 4
                        * gamma
                        * LSinv ^ 2
                        * pole_disp
                        * r ^ 2
                        * tc
                        * (
                            1
                            + 6 * gamma * pole_disp * tc
                            + 6 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                        + 2
                        * r
                        * (LSinv - LTinv + dLSinvdr * r - dLTinvdr * r)
                        * (
                            1
                            + 5 * gamma * pole_disp * tc
                            + 9 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                            + 6 * gamma ^ 3 * pole_disp ^ 3 * tc ^ 3
                        )
                    )
                )
                + D
                * (
                    -(
                        cs ^ 2
                        * LSinv
                        * r ^ 2
                        * tc
                        * (
                            dDdr
                            + 2 * (2 * dDdr + dNr2dr) * gamma * pole_disp * tc
                            + 10
                            * (dDdr + dNr2dr)
                            * gamma ^ 2
                            * pole_disp ^ 2
                            * tc ^ 2
                            + 6
                            * (2 * dDdr + 3 * dNr2dr)
                            * gamma ^ 3
                            * pole_disp ^ 3
                            * tc ^ 3
                            + 6
                            * (dDdr + 2 * dNr2dr)
                            * gamma ^ 4
                            * pole_disp ^ 4
                            * tc ^ 4
                        )
                    )
                    - D ^ 2
                    * tc
                    * (
                        (-1 + gamma) * r ^ 2
                        + 4 * (-1 + gamma) * gamma * pole_disp * r ^ 2 * tc
                        + 10
                        * (-1 + gamma)
                        * gamma ^ 2
                        * pole_disp ^ 2
                        * r ^ 2
                        * tc ^ 2
                        + 6
                        * gamma ^ 3
                        * pole_disp
                        * (
                            cs ^ 2 * m ^ 2
                            + 2 * (-1 + gamma) * pole_disp ^ 2 * r ^ 2
                        )
                        * tc ^ 3
                        + 6
                        * gamma ^ 4
                        * pole_disp ^ 2
                        * (
                            2 * cs ^ 2 * m ^ 2
                            + (-1 + gamma) * pole_disp ^ 2 * r ^ 2
                        )
                        * tc ^ 4
                        + 6 * cs ^ 2 * gamma ^ 5 * m ^ 2 * pole_disp ^ 3 * tc ^ 5
                    )
                    + cs ^ 2
                    * D
                    * (
                        2
                        * gamma
                        * LSinv ^ 2
                        * pole_disp
                        * r ^ 2
                        * tc ^ 2
                        * (
                            1
                            + 5 * gamma * pole_disp * tc
                            + 9 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                            + 6 * gamma ^ 3 * pole_disp ^ 3 * tc ^ 3
                        )
                        - dtcdr
                        * LSinv
                        * r ^ 2
                        * (
                            -1
                            - 2 * gamma * pole_disp * tc
                            + 6 * gamma ^ 3 * pole_disp ^ 3 * tc ^ 3
                            + 6 * gamma ^ 4 * pole_disp ^ 4 * tc ^ 4
                        )
                        + r
                        * (LSinv - LTinv + dLSinvdr * r - dLTinvdr * r)
                        * tc
                        * (
                            1
                            + 4 * gamma * pole_disp * tc
                            + 10 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                            + 12 * gamma ^ 3 * pole_disp ^ 3 * tc ^ 3
                            + 6 * gamma ^ 4 * pole_disp ^ 4 * tc ^ 4
                        )
                    )
                )
            )
            + gamma
            * omega_tilde ^ 3
            * (
                -(
                    gamma ^ 2
                    * Nr2 ^ 3
                    * pole_disp ^ 2
                    * tc ^ 3
                    * (
                        gamma
                        * pole_disp ^ 2
                        * r ^ 2
                        * (
                            -3
                            - 6 * pole_disp * tc
                            + 4 * (-2 + gamma) * gamma * pole_disp ^ 2 * tc ^ 2
                        )
                        + cs ^ 2
                        * m ^ 2
                        * (
                            2
                            + 8 * gamma * pole_disp * tc
                            + 11 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                            + 6 * gamma ^ 3 * pole_disp ^ 3 * tc ^ 3
                        )
                    )
                )
                - 2
                * D
                * pole_disp
                * (1 + gamma * pole_disp * tc)
                * (
                    D ^ 2
                    * tc
                    * (1 + gamma * pole_disp * tc)
                    * (
                        (-1 + gamma) * pole_disp * r ^ 2
                        + gamma
                        * (
                            3 * cs ^ 2 * m ^ 2
                            + 2 * (-1 + gamma) * pole_disp ^ 2 * r ^ 2
                        )
                        * tc
                        + 2
                        * gamma ^ 2
                        * pole_disp
                        * (
                            3 * cs ^ 2 * m ^ 2
                            + (-1 + gamma) * pole_disp ^ 2 * r ^ 2
                        )
                        * tc ^ 2
                        + 3 * cs ^ 2 * gamma ^ 3 * m ^ 2 * pole_disp ^ 2 * tc ^ 3
                    )
                    + cs ^ 2
                    * LSinv
                    * pole_disp
                    * r ^ 2
                    * tc
                    * (
                        dDdr
                        + (3 * dDdr + 2 * dNr2dr) * gamma * pole_disp * tc
                        + (4 * dDdr + 5 * dNr2dr)
                        * gamma ^ 2
                        * pole_disp ^ 2
                        * tc ^ 2
                        + 2
                        * (dDdr + 2 * dNr2dr)
                        * gamma ^ 3
                        * pole_disp ^ 3
                        * tc ^ 3
                    )
                    - cs ^ 2
                    * D
                    * pole_disp
                    * (
                        gamma
                        * LSinv ^ 2
                        * pole_disp
                        * r ^ 2
                        * tc ^ 2
                        * (
                            2
                            + 5 * gamma * pole_disp * tc
                            + 4 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                        - dtcdr
                        * LSinv
                        * r ^ 2
                        * (
                            -1
                            - gamma * pole_disp * tc
                            + gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                            + 2 * gamma ^ 3 * pole_disp ^ 3 * tc ^ 3
                        )
                        + r
                        * (LSinv - LTinv + dLSinvdr * r - dLTinvdr * r)
                        * tc
                        * (
                            1
                            + 3 * gamma * pole_disp * tc
                            + 4 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                            + 2 * gamma ^ 3 * pole_disp ^ 3 * tc ^ 3
                        )
                    )
                )
                + gamma
                * Nr2 ^ 2
                * pole_disp
                * tc ^ 2
                * (
                    cs ^ 2
                    * gamma
                    * pole_disp ^ 3
                    * tc
                    * (
                        2
                        * gamma
                        * LSinv ^ 2
                        * pole_disp
                        * r ^ 2
                        * tc
                        * (3 + 4 * gamma * pole_disp * tc)
                        + r
                        * (LSinv - LTinv + dLSinvdr * r - dLTinvdr * r)
                        * (
                            3
                            + 6 * gamma * pole_disp * tc
                            + 4 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                    )
                    - D
                    * (
                        gamma
                        * pole_disp ^ 2
                        * r ^ 2
                        * (
                            -4
                            - (11 + 3 * gamma) * pole_disp * tc
                            + 6 * gamma * (-5 + 2 * gamma) * pole_disp ^ 2 * tc ^ 2
                            + 4
                            * gamma ^ 2
                            * (-5 + 3 * gamma)
                            * pole_disp ^ 3
                            * tc ^ 3
                        )
                        + 2
                        * cs ^ 2
                        * m ^ 2
                        * (
                            2
                            + 10 * gamma * pole_disp * tc
                            + 22 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                            + 23 * gamma ^ 3 * pole_disp ^ 3 * tc ^ 3
                            + 9 * gamma ^ 4 * pole_disp ^ 4 * tc ^ 4
                        )
                    )
                )
                + Nr2
                * tc
                * (
                    -(
                        cs ^ 2
                        * dNr2dr
                        * gamma ^ 2
                        * LSinv
                        * pole_disp ^ 4
                        * r ^ 2
                        * tc ^ 2
                        * (
                            3
                            + 6 * gamma * pole_disp * tc
                            + 4 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                    )
                    + cs ^ 2
                    * D
                    * gamma
                    * pole_disp ^ 3
                    * tc
                    * (
                        8
                        * gamma
                        * LSinv ^ 2
                        * pole_disp
                        * r ^ 2
                        * tc
                        * (
                            1
                            + 3 * gamma * pole_disp * tc
                            + 2 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                        - dtcdr
                        * gamma
                        * LSinv
                        * pole_disp
                        * r ^ 2
                        * (
                            3
                            + 6 * gamma * pole_disp * tc
                            + 4 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                        + 2
                        * r
                        * (LSinv - LTinv + dLSinvdr * r - dLTinvdr * r)
                        * (
                            2
                            + 7 * gamma * pole_disp * tc
                            + 9 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                            + 4 * gamma ^ 3 * pole_disp ^ 3 * tc ^ 3
                        )
                    )
                    - D ^ 2
                    * (1 + gamma * pole_disp * tc)
                    * (
                        2
                        * gamma
                        * pole_disp ^ 2
                        * r ^ 2
                        * (
                            -1
                            + (-4 + gamma) * pole_disp * tc
                            + 2 * gamma * (-5 + 3 * gamma) * pole_disp ^ 2 * tc ^ 2
                            + 2
                            * gamma ^ 2
                            * (-4 + 3 * gamma)
                            * pole_disp ^ 3
                            * tc ^ 3
                        )
                        + cs ^ 2
                        * m ^ 2
                        * (
                            1
                            + 9 * gamma * pole_disp * tc
                            + 31 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                            + 41 * gamma ^ 3 * pole_disp ^ 3 * tc ^ 3
                            + 18 * gamma ^ 4 * pole_disp ^ 4 * tc ^ 4
                        )
                    )
                )
            )
            - omega_tilde
            * pole_disp
            * (1 + gamma * pole_disp * tc)
            * (
                gamma
                * pole_disp ^ 3
                * r ^ 2
                * tc
                * (gamma * Nr2 * pole_disp * tc + D * (1 + gamma * pole_disp * tc)) ^ 2
                * (
                    gamma * Nr2 * (-1 + (-2 + gamma) * pole_disp * tc)
                    + D * (-1 + gamma) * (1 + gamma * pole_disp * tc)
                )
                + cs ^ 2
                * (
                    2 * D ^ 3 * m ^ 2 * (1 + gamma * pole_disp * tc) ^ 5
                    + gamma ^ 3
                    * Nr2
                    * pole_disp ^ 3
                    * tc ^ 3
                    * (
                        dNr2dr
                        * LSinv
                        * pole_disp ^ 2
                        * r ^ 2
                        * (1 + gamma * pole_disp * tc)
                        + m ^ 2
                        * Nr2 ^ 2
                        * (
                            1
                            + 3 * gamma * pole_disp * tc
                            + 2 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                        + Nr2
                        * pole_disp ^ 2
                        * (
                            -2 * gamma * LSinv ^ 2 * pole_disp * r ^ 2 * tc
                            - r
                            * (LSinv - LTinv + dLSinvdr * r - dLTinvdr * r)
                            * (1 + gamma * pole_disp * tc)
                        )
                    )
                    + D
                    * gamma
                    * pole_disp ^ 2
                    * tc
                    * (1 + gamma * pole_disp * tc)
                    * (
                        LSinv
                        * pole_disp
                        * r ^ 2
                        * (1 + gamma * pole_disp * tc)
                        * (dDdr + (dDdr + 2 * dNr2dr) * gamma * pole_disp * tc)
                        + 2
                        * gamma
                        * m ^ 2
                        * Nr2 ^ 2
                        * tc
                        * (
                            2
                            + 5 * gamma * pole_disp * tc
                            + 3 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                        - gamma
                        * Nr2
                        * pole_disp ^ 2
                        * tc
                        * (
                            -(dtcdr * gamma * LSinv * pole_disp * r ^ 2)
                            + 4 * gamma * LSinv ^ 2 * pole_disp * r ^ 2 * tc
                            + 2
                            * r
                            * (LSinv - LTinv + dLSinvdr * r - dLTinvdr * r)
                            * (1 + gamma * pole_disp * tc)
                        )
                    )
                    + D ^ 2
                    * gamma
                    * pole_disp
                    * (1 + gamma * pole_disp * tc) ^ 2
                    * (
                        m ^ 2
                        * Nr2
                        * tc
                        * (
                            5
                            + 11 * gamma * pole_disp * tc
                            + 6 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                        )
                        + pole_disp ^ 2
                        * (
                            -2 * gamma * LSinv ^ 2 * pole_disp * r ^ 2 * tc ^ 2
                            - dtcdr * LSinv * r ^ 2 * (1 - gamma * pole_disp * tc)
                            - r
                            * (LSinv - LTinv + dLSinvdr * r - dLTinvdr * r)
                            * tc
                            * (1 + gamma * pole_disp * tc)
                        )
                    )
                )
            )
        )
        + cs ^ 2
        * LTinv
        * (omega_tilde ^ 2 + pole_disp ^ 2)
        * r
        * (
            -(
                dSigmadr
                * gamma
                * omega_tilde
                * (omega_tilde ^ 2 + pole_disp ^ 2)
                * r
                * tc
                * (
                    1
                    + 2 * gamma * pole_disp * tc
                    + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                )
                * (
                    gamma ^ 2
                    * Nr2 ^ 2
                    * (omega_tilde ^ 2 + pole_disp ^ 2)
                    * tc ^ 2
                    + 2
                    * D
                    * gamma
                    * Nr2
                    * tc
                    * (pole_disp + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc)
                    + D ^ 2
                    * (
                        1
                        + 2 * gamma * pole_disp * tc
                        + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                    )
                )
            )
            + Sigma
            * (
                2
                * m
                * Omega
                * (
                    pole_disp
                    + gamma * (-(omega_tilde ^ 2) + 3 * pole_disp ^ 2) * tc
                    + gamma ^ 2
                    * pole_disp
                    * (-(omega_tilde ^ 2) + 3 * pole_disp ^ 2)
                    * tc ^ 2
                    + gamma ^ 3 * (-(omega_tilde ^ 4) + pole_disp ^ 4) * tc ^ 3
                )
                * (
                    gamma ^ 2
                    * Nr2 ^ 2
                    * (omega_tilde ^ 2 + pole_disp ^ 2)
                    * tc ^ 2
                    + 2
                    * D
                    * gamma
                    * Nr2
                    * tc
                    * (pole_disp + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc)
                    + D ^ 2
                    * (
                        1
                        + 2 * gamma * pole_disp * tc
                        + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                    )
                )
                + gamma
                * (omega_tilde ^ 2 + pole_disp ^ 2)
                * (
                    -(
                        LSinv
                        * omega_tilde
                        * r
                        * tc
                        * (
                            -1
                            + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                        )
                        * (
                            gamma ^ 2
                            * Nr2 ^ 2
                            * (omega_tilde ^ 2 + pole_disp ^ 2)
                            * tc ^ 2
                            + 2
                            * D
                            * gamma
                            * Nr2
                            * tc
                            * (
                                pole_disp
                                + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc
                            )
                            + D ^ 2
                            * (
                                1
                                + 2 * gamma * pole_disp * tc
                                + gamma ^ 2
                                * (omega_tilde ^ 2 + pole_disp ^ 2)
                                * tc ^ 2
                            )
                        )
                    )
                    + r
                    * (
                        1
                        + 2 * gamma * pole_disp * tc
                        + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                    )
                    * (
                        -(
                            dOmegadr
                            * gamma ^ 2
                            * m
                            * (D + Nr2) ^ 2
                            * omega_tilde ^ 2
                            * tc ^ 3
                        )
                        + gamma ^ 2
                        * (D + Nr2)
                        * omega_tilde ^ 3
                        * tc ^ 2
                        * (dtcdr * (D + Nr2) + (dDdr + dNr2dr) * tc)
                        + dOmegadr
                        * m
                        * (D + Nr2)
                        * tc
                        * (1 + gamma * pole_disp * tc)
                        * (
                            gamma * Nr2 * pole_disp * tc
                            + D * (1 + gamma * pole_disp * tc)
                        )
                        + omega_tilde
                        * (
                            D ^ 2
                            * dtcdr
                            * (-1 + gamma ^ 2 * pole_disp ^ 2 * tc ^ 2)
                            + Nr2
                            * tc
                            * (
                                dDdr
                                + gamma
                                * pole_disp
                                * (2 * dDdr + dtcdr * gamma * Nr2 * pole_disp)
                                * tc
                                + (dDdr + dNr2dr)
                                * gamma ^ 2
                                * pole_disp ^ 2
                                * tc ^ 2
                            )
                            + D
                            * (
                                -(dtcdr * Nr2)
                                + (dDdr - dNr2dr) * tc
                                + 2
                                * gamma
                                * pole_disp
                                * (dDdr + dtcdr * gamma * Nr2 * pole_disp)
                                * tc ^ 2
                                + (dDdr + dNr2dr)
                                * gamma ^ 2
                                * pole_disp ^ 2
                                * tc ^ 3
                            )
                        )
                    )
                )
            )
        )
    ) / (
        cs ^ 2
        * (omega_tilde ^ 2 + pole_disp ^ 2) ^ 2
        * r ^ 2
        * Sigma
        * (
            1
            + 2 * gamma * pole_disp * tc
            + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
        )
        ^ 2
        * (
            gamma ^ 2 * Nr2 ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
            + 2
            * D
            * gamma
            * Nr2
            * tc
            * (pole_disp + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc)
            + D ^ 2
            * (
                1
                + 2 * gamma * pole_disp * tc
                + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
            )
        )
    )

    PsiR = (
        -(
            (
                dSigmadr
                * (omega_tilde ^ 2 + pole_disp ^ 2)
                * r
                * (
                    -2 * m * Omega * omega_tilde * Phi
                    + dPhidr * (omega_tilde ^ 2 + pole_disp ^ 2) * r
                )
            )
            / Sigma
        )
        + (
            (omega_tilde ^ 2 + pole_disp ^ 2) ^ 2
            * r
            * (
                D ^ 2
                * (
                    1
                    + 2 * gamma * pole_disp * tc
                    + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                )
                * (
                    -(
                        d2Phidr2
                        * r
                        * (
                            1
                            + 2 * gamma * pole_disp * tc
                            + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                        )
                    )
                    + dPhidr
                    * (
                        -1
                        + gamma * pole_disp * (-2 - LSinv * r) * tc
                        + gamma ^ 2
                        * (omega_tilde ^ 2 + pole_disp ^ 2)
                        * (-1 - LSinv * r)
                        * tc ^ 2
                    )
                )
                + gamma
                * Nr2
                * tc
                * (
                    dPhidr
                    * r
                    * (
                        dDdr * pole_disp
                        + (dDdr + dNr2dr)
                        * gamma
                        * (omega_tilde ^ 2 + pole_disp ^ 2)
                        * tc
                    )
                    * (
                        1
                        + 2 * gamma * pole_disp * tc
                        + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                    )
                    + gamma
                    * Nr2
                    * (
                        dPhidr * dtcdr * (omega_tilde ^ 2 + pole_disp ^ 2) * r
                        + (
                            -(d2Phidr2 * (omega_tilde ^ 2 + pole_disp ^ 2) * r)
                            + dPhidr
                            * (
                                -(dOmegadr * m * omega_tilde * r)
                                + omega_tilde ^ 2
                                * (-1 + dtcdr * gamma * pole_disp * r)
                                + pole_disp ^ 2 * (-1 + dtcdr * gamma * pole_disp * r)
                            )
                        )
                        * tc
                        + gamma
                        * pole_disp
                        * (
                            -2 * d2Phidr2 * (omega_tilde ^ 2 + pole_disp ^ 2) * r
                            + dPhidr
                            * (
                                -(LSinv * (omega_tilde ^ 2 + pole_disp ^ 2) * r)
                                - 2
                                * (
                                    omega_tilde ^ 2
                                    + pole_disp ^ 2
                                    + dOmegadr * m * omega_tilde * r
                                )
                            )
                        )
                        * tc ^ 2
                        + gamma ^ 2
                        * (omega_tilde ^ 2 + pole_disp ^ 2) ^ 2
                        * (-(d2Phidr2 * r) + dPhidr * (-1 - LSinv * r))
                        * tc ^ 3
                    )
                )
                + D
                * (
                    dPhidr
                    * r
                    * (
                        1
                        + 2 * gamma * pole_disp * tc
                        + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                    )
                    * (
                        dDdr
                        + (2 * dDdr + dNr2dr) * gamma * pole_disp * tc
                        + (dDdr + dNr2dr)
                        * gamma ^ 2
                        * (omega_tilde ^ 2 + pole_disp ^ 2)
                        * tc ^ 2
                    )
                    + gamma
                    * Nr2
                    * (
                        dPhidr * dtcdr * pole_disp * r
                        + 2
                        * (
                            -(d2Phidr2 * pole_disp * r)
                            + dPhidr
                            * (
                                -pole_disp
                                + dtcdr
                                * gamma
                                * (omega_tilde ^ 2 + pole_disp ^ 2)
                                * r
                            )
                        )
                        * tc
                        + gamma
                        * (
                            -2 * d2Phidr2 * (omega_tilde ^ 2 + 3 * pole_disp ^ 2) * r
                            + dPhidr
                            * (
                                -2 * dOmegadr * m * omega_tilde * r
                                - 2 * LSinv * pole_disp ^ 2 * r
                                + pole_disp ^ 2 * (-6 + dtcdr * gamma * pole_disp * r)
                                + omega_tilde ^ 2
                                * (-2 + dtcdr * gamma * pole_disp * r)
                            )
                        )
                        * tc ^ 2
                        - 2
                        * gamma ^ 2
                        * pole_disp
                        * (
                            3 * d2Phidr2 * (omega_tilde ^ 2 + pole_disp ^ 2) * r
                            + dPhidr
                            * (
                                3 * omega_tilde ^ 2
                                + 3 * pole_disp ^ 2
                                + dOmegadr * m * omega_tilde * r
                                + 2 * LSinv * (omega_tilde ^ 2 + pole_disp ^ 2) * r
                            )
                        )
                        * tc ^ 3
                        + 2
                        * gamma ^ 3
                        * (omega_tilde ^ 2 + pole_disp ^ 2) ^ 2
                        * (-(d2Phidr2 * r) + dPhidr * (-1 - LSinv * r))
                        * tc ^ 4
                    )
                )
            )
            - m
            * Phi
            * (
                4
                * m
                * Omega ^ 2
                * (-(omega_tilde ^ 2) + pole_disp ^ 2)
                * (
                    1
                    + 2 * gamma * pole_disp * tc
                    + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                )
                * (
                    gamma ^ 2
                    * Nr2 ^ 2
                    * (omega_tilde ^ 2 + pole_disp ^ 2)
                    * tc ^ 2
                    + 2
                    * D
                    * gamma
                    * Nr2
                    * tc
                    * (pole_disp + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc)
                    + D ^ 2
                    * (
                        1
                        + 2 * gamma * pole_disp * tc
                        + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                    )
                )
                - (
                    gamma ^ 2
                    * Nr2 ^ 2
                    * (omega_tilde ^ 2 + pole_disp ^ 2)
                    * tc ^ 2
                    + 2
                    * D
                    * gamma
                    * Nr2
                    * tc
                    * (pole_disp + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc)
                    + D ^ 2
                    * (
                        1
                        + 2 * gamma * pole_disp * tc
                        + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                    )
                )
                * (
                    D
                    * m
                    * (-(omega_tilde ^ 2) + pole_disp ^ 2)
                    * (
                        1
                        + 2 * gamma * pole_disp * tc
                        + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                    )
                    + (omega_tilde ^ 2 + pole_disp ^ 2)
                    * (
                        gamma
                        * m
                        * Nr2
                        * tc
                        * (
                            pole_disp
                            + gamma * (-(omega_tilde ^ 2) + pole_disp ^ 2) * tc
                        )
                        + 2
                        * dOmegadr
                        * omega_tilde
                        * r
                        * (
                            1
                            + 2 * gamma * pole_disp * tc
                            + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                        )
                    )
                )
                + 2
                * Omega
                * (
                    -(
                        gamma ^ 2
                        * LSinv
                        * omega_tilde
                        * (omega_tilde ^ 2 + pole_disp ^ 2) ^ 2
                        * r
                        * tc ^ 2
                        * (
                            gamma ^ 2
                            * Nr2 ^ 2
                            * (omega_tilde ^ 2 + pole_disp ^ 2)
                            * tc ^ 2
                            + 2
                            * D
                            * gamma
                            * Nr2
                            * tc
                            * (
                                pole_disp
                                + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc
                            )
                            + D ^ 2
                            * (
                                1
                                + 2 * gamma * pole_disp * tc
                                + gamma ^ 2
                                * (omega_tilde ^ 2 + pole_disp ^ 2)
                                * tc ^ 2
                            )
                        )
                    )
                    + r
                    * (
                        gamma ^ 2
                        * Nr2 ^ 2
                        * (omega_tilde ^ 2 + pole_disp ^ 2)
                        * tc
                        * (
                            dtcdr * omega_tilde * (omega_tilde ^ 2 + pole_disp ^ 2)
                            + 2
                            * (
                                dOmegadr * m * (-(omega_tilde ^ 2) + pole_disp ^ 2)
                                + dtcdr
                                * gamma
                                * omega_tilde
                                * pole_disp
                                * (omega_tilde ^ 2 + pole_disp ^ 2)
                            )
                            * tc
                            + dOmegadr
                            * gamma
                            * m
                            * pole_disp
                            * (-5 * omega_tilde ^ 2 + 3 * pole_disp ^ 2)
                            * tc ^ 2
                            + dOmegadr
                            * gamma ^ 2
                            * m
                            * (-(omega_tilde ^ 4) + pole_disp ^ 4)
                            * tc ^ 3
                        )
                        + D
                        * (
                            1
                            + 2 * gamma * pole_disp * tc
                            + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                        )
                        * (
                            -(
                                D
                                * dOmegadr
                                * gamma ^ 2
                                * m
                                * omega_tilde ^ 4
                                * tc ^ 2
                            )
                            + (dDdr + dNr2dr) * gamma ^ 2 * omega_tilde ^ 5 * tc ^ 2
                            + D
                            * dOmegadr
                            * m
                            * pole_disp ^ 2
                            * (1 + gamma * pole_disp * tc) ^ 2
                            - D
                            * dOmegadr
                            * m
                            * omega_tilde ^ 2
                            * (1 + 2 * gamma * pole_disp * tc)
                            + omega_tilde
                            * pole_disp ^ 2
                            * (
                                dDdr
                                + 2 * dDdr * gamma * pole_disp * tc
                                + (dDdr + dNr2dr)
                                * gamma ^ 2
                                * pole_disp ^ 2
                                * tc ^ 2
                            )
                            + omega_tilde ^ 3
                            * (
                                dDdr
                                + 2 * dDdr * gamma * pole_disp * tc
                                + 2
                                * (dDdr + dNr2dr)
                                * gamma ^ 2
                                * pole_disp ^ 2
                                * tc ^ 2
                            )
                        )
                        + gamma
                        * Nr2
                        * tc
                        * (
                            -2
                            * D
                            * dOmegadr
                            * gamma ^ 3
                            * m
                            * omega_tilde ^ 6
                            * tc ^ 3
                            + (dDdr + dNr2dr) * gamma ^ 3 * omega_tilde ^ 7 * tc ^ 3
                            + D
                            * dOmegadr
                            * m
                            * pole_disp ^ 3
                            * (1 + gamma * pole_disp * tc) ^ 2
                            * (3 + 2 * gamma * pole_disp * tc)
                            - D
                            * dOmegadr
                            * gamma
                            * m
                            * omega_tilde ^ 4
                            * tc
                            * (
                                4
                                + 9 * gamma * pole_disp * tc
                                + 2 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                            )
                            + omega_tilde
                            * pole_disp ^ 3
                            * (1 + gamma * pole_disp * tc)
                            * (
                                2 * (dDdr + D * dtcdr * gamma * pole_disp)
                                + (3 * dDdr + dNr2dr) * gamma * pole_disp * tc
                                + (dDdr + dNr2dr)
                                * gamma ^ 2
                                * pole_disp ^ 2
                                * tc ^ 2
                            )
                            + gamma
                            * omega_tilde ^ 5
                            * (
                                2 * D * dtcdr
                                + (dDdr + dNr2dr + 2 * D * dtcdr * gamma * pole_disp)
                                * tc
                                + 2 * (2 * dDdr + dNr2dr) * gamma * pole_disp * tc ^ 2
                                + 3
                                * (dDdr + dNr2dr)
                                * gamma ^ 2
                                * pole_disp ^ 2
                                * tc ^ 3
                            )
                            + D
                            * dOmegadr
                            * m
                            * omega_tilde ^ 2
                            * pole_disp
                            * (
                                -1
                                - 4 * gamma * pole_disp * tc
                                - 2 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                                + 2 * gamma ^ 3 * pole_disp ^ 3 * tc ^ 3
                            )
                            + omega_tilde ^ 3
                            * pole_disp
                            * (
                                2 * (dDdr + 2 * D * dtcdr * gamma * pole_disp)
                                + 2
                                * gamma
                                * pole_disp
                                * (
                                    3 * dDdr
                                    + dNr2dr
                                    + 2 * D * dtcdr * gamma * pole_disp
                                )
                                * tc
                                + 4
                                * (2 * dDdr + dNr2dr)
                                * gamma ^ 2
                                * pole_disp ^ 2
                                * tc ^ 2
                                + 3
                                * (dDdr + dNr2dr)
                                * gamma ^ 3
                                * pole_disp ^ 3
                                * tc ^ 3
                            )
                        )
                    )
                )
            )
        )
        / (
            (
                1
                + 2 * gamma * pole_disp * tc
                + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
            )
            * (
                gamma ^ 2 * Nr2 ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                + 2
                * D
                * gamma
                * Nr2
                * tc
                * (pole_disp + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc)
                + D ^ 2
                * (
                    1
                    + 2 * gamma * pole_disp * tc
                    + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                )
            )
        )
    ) / ((omega_tilde ^ 2 + pole_disp ^ 2) ^ 2 * r ^ 2)

    PsiI = (
        (
            -2
            * dSigmadr
            * m
            * Omega
            * Phi
            * pole_disp
            * (omega_tilde ^ 2 + pole_disp ^ 2)
            * r
        )
        / Sigma
        + (
            m
            * Phi
            * (
                -8
                * m
                * Omega ^ 2
                * omega_tilde
                * pole_disp
                * (
                    1
                    + 2 * gamma * pole_disp * tc
                    + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                )
                * (
                    gamma ^ 2
                    * Nr2 ^ 2
                    * (omega_tilde ^ 2 + pole_disp ^ 2)
                    * tc ^ 2
                    + 2
                    * D
                    * gamma
                    * Nr2
                    * tc
                    * (pole_disp + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc)
                    + D ^ 2
                    * (
                        1
                        + 2 * gamma * pole_disp * tc
                        + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                    )
                )
                + (
                    gamma ^ 2
                    * Nr2 ^ 2
                    * (omega_tilde ^ 2 + pole_disp ^ 2)
                    * tc ^ 2
                    + 2
                    * D
                    * gamma
                    * Nr2
                    * tc
                    * (pole_disp + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc)
                    + D ^ 2
                    * (
                        1
                        + 2 * gamma * pole_disp * tc
                        + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                    )
                )
                * (
                    -2
                    * dOmegadr
                    * gamma ^ 2
                    * omega_tilde ^ 4
                    * pole_disp
                    * r
                    * tc ^ 2
                    - 2
                    * dOmegadr
                    * pole_disp ^ 3
                    * r
                    * (1 + gamma * pole_disp * tc) ^ 2
                    + gamma
                    * m
                    * omega_tilde ^ 3
                    * tc
                    * (
                        Nr2
                        + 2 * D * gamma * pole_disp * tc
                        + 2 * gamma * Nr2 * pole_disp * tc
                    )
                    - 2
                    * dOmegadr
                    * omega_tilde ^ 2
                    * pole_disp
                    * r
                    * (
                        1
                        + 2 * gamma * pole_disp * tc
                        + 2 * gamma ^ 2 * pole_disp ^ 2 * tc ^ 2
                    )
                    + m
                    * omega_tilde
                    * pole_disp
                    * (
                        2 * D * (1 + gamma * pole_disp * tc) ^ 2
                        + gamma
                        * Nr2
                        * pole_disp
                        * tc
                        * (1 + 2 * gamma * pole_disp * tc)
                    )
                )
                + 2
                * Omega
                * (
                    -(
                        gamma
                        * LSinv
                        * (omega_tilde ^ 2 + pole_disp ^ 2) ^ 2
                        * r
                        * tc
                        * (1 + gamma * pole_disp * tc)
                        * (
                            gamma ^ 2
                            * Nr2 ^ 2
                            * (omega_tilde ^ 2 + pole_disp ^ 2)
                            * tc ^ 2
                            + 2
                            * D
                            * gamma
                            * Nr2
                            * tc
                            * (
                                pole_disp
                                + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc
                            )
                            + D ^ 2
                            * (
                                1
                                + 2 * gamma * pole_disp * tc
                                + gamma ^ 2
                                * (omega_tilde ^ 2 + pole_disp ^ 2)
                                * tc ^ 2
                            )
                        )
                    )
                    + r
                    * (
                        -2
                        * D ^ 2
                        * dOmegadr
                        * m
                        * omega_tilde
                        * pole_disp
                        * (
                            1
                            + 2 * gamma * pole_disp * tc
                            + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                        )
                        ^ 2
                        + gamma
                        * Nr2
                        * (omega_tilde ^ 2 + pole_disp ^ 2)
                        * tc
                        * (
                            -(dDdr * omega_tilde ^ 2)
                            + dtcdr * gamma * Nr2 * omega_tilde ^ 2 * pole_disp
                            + dDdr * pole_disp ^ 2
                            + dtcdr * gamma * Nr2 * pole_disp ^ 3
                            + gamma
                            * (
                                dNr2dr * pole_disp * (omega_tilde ^ 2 + pole_disp ^ 2)
                                + dDdr
                                * pole_disp
                                * (-(omega_tilde ^ 2) + 3 * pole_disp ^ 2)
                                + Nr2
                                * (
                                    -(dtcdr * gamma * omega_tilde ^ 4)
                                    - 4 * dOmegadr * m * omega_tilde * pole_disp
                                    + dtcdr * gamma * pole_disp ^ 4
                                )
                            )
                            * tc
                            + gamma ^ 2
                            * (
                                dOmegadr
                                * m
                                * Nr2
                                * omega_tilde
                                * (omega_tilde ^ 2 - 7 * pole_disp ^ 2)
                                + 2
                                * dNr2dr
                                * pole_disp ^ 2
                                * (omega_tilde ^ 2 + pole_disp ^ 2)
                                + dDdr
                                * (
                                    -(omega_tilde ^ 4)
                                    + 2 * omega_tilde ^ 2 * pole_disp ^ 2
                                    + 3 * pole_disp ^ 4
                                )
                            )
                            * tc ^ 2
                            + gamma ^ 3
                            * pole_disp
                            * (omega_tilde ^ 2 + pole_disp ^ 2)
                            * (
                                -2 * dOmegadr * m * Nr2 * omega_tilde
                                + dDdr * (omega_tilde ^ 2 + pole_disp ^ 2)
                                + dNr2dr * (omega_tilde ^ 2 + pole_disp ^ 2)
                            )
                            * tc ^ 3
                        )
                        + D
                        * (
                            (omega_tilde ^ 2 + pole_disp ^ 2)
                            * (
                                dDdr * pole_disp
                                + dtcdr
                                * gamma
                                * Nr2
                                * (omega_tilde ^ 2 + pole_disp ^ 2)
                            )
                            + gamma
                            * (
                                4
                                * dDdr
                                * pole_disp ^ 2
                                * (omega_tilde ^ 2 + pole_disp ^ 2)
                                + dNr2dr * (omega_tilde ^ 2 + pole_disp ^ 2) ^ 2
                                + Nr2
                                * (
                                    -(dOmegadr * m * omega_tilde ^ 3)
                                    + 2 * dtcdr * gamma * omega_tilde ^ 4 * pole_disp
                                    - 5 * dOmegadr * m * omega_tilde * pole_disp ^ 2
                                    + 4
                                    * dtcdr
                                    * gamma
                                    * omega_tilde ^ 2
                                    * pole_disp ^ 3
                                    + 2 * dtcdr * gamma * pole_disp ^ 5
                                )
                            )
                            * tc
                            + gamma ^ 2
                            * (
                                3
                                * dNr2dr
                                * pole_disp
                                * (omega_tilde ^ 2 + pole_disp ^ 2) ^ 2
                                + 2
                                * dDdr
                                * pole_disp
                                * (
                                    omega_tilde ^ 4
                                    + 4 * omega_tilde ^ 2 * pole_disp ^ 2
                                    + 3 * pole_disp ^ 4
                                )
                                - Nr2
                                * (
                                    dtcdr * gamma * omega_tilde ^ 6
                                    + 8 * dOmegadr * m * omega_tilde ^ 3 * pole_disp
                                    + dtcdr * gamma * omega_tilde ^ 4 * pole_disp ^ 2
                                    + 16 * dOmegadr * m * omega_tilde * pole_disp ^ 3
                                    - dtcdr * gamma * omega_tilde ^ 2 * pole_disp ^ 4
                                    - dtcdr * gamma * pole_disp ^ 6
                                )
                            )
                            * tc ^ 2
                            + gamma ^ 3
                            * (omega_tilde ^ 2 + pole_disp ^ 2)
                            * (
                                dOmegadr
                                * m
                                * Nr2
                                * omega_tilde
                                * (omega_tilde ^ 2 - 15 * pole_disp ^ 2)
                                + 4
                                * dDdr
                                * pole_disp ^ 2
                                * (omega_tilde ^ 2 + pole_disp ^ 2)
                                + dNr2dr
                                * (
                                    omega_tilde ^ 4
                                    + 4 * omega_tilde ^ 2 * pole_disp ^ 2
                                    + 3 * pole_disp ^ 4
                                )
                            )
                            * tc ^ 3
                            + gamma ^ 4
                            * pole_disp
                            * (omega_tilde ^ 2 + pole_disp ^ 2) ^ 2
                            * (
                                -4 * dOmegadr * m * Nr2 * omega_tilde
                                + dDdr * (omega_tilde ^ 2 + pole_disp ^ 2)
                                + dNr2dr * (omega_tilde ^ 2 + pole_disp ^ 2)
                            )
                            * tc ^ 4
                        )
                    )
                )
            )
            - dPhidr
            * gamma
            * (omega_tilde ^ 2 + pole_disp ^ 2) ^ 2
            * r
            * (
                -(
                    LSinv
                    * omega_tilde
                    * r
                    * tc
                    * (
                        gamma ^ 2
                        * Nr2 ^ 2
                        * (omega_tilde ^ 2 + pole_disp ^ 2)
                        * tc ^ 2
                        + 2
                        * D
                        * gamma
                        * Nr2
                        * tc
                        * (pole_disp + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc)
                        + D ^ 2
                        * (
                            1
                            + 2 * gamma * pole_disp * tc
                            + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                        )
                    )
                )
                + r
                * (
                    dOmegadr
                    * gamma ^ 2
                    * m
                    * Nr2
                    * (D + Nr2)
                    * omega_tilde ^ 2
                    * tc ^ 3
                    + gamma ^ 2
                    * omega_tilde ^ 3
                    * tc ^ 2
                    * (
                        -(Nr2 * (dtcdr * Nr2 + dDdr * tc))
                        + D * (-(dtcdr * Nr2) + dNr2dr * tc)
                    )
                    - dOmegadr
                    * m
                    * Nr2
                    * tc
                    * (1 + gamma * pole_disp * tc)
                    * (gamma * Nr2 * pole_disp * tc + D * (1 + gamma * pole_disp * tc))
                    + omega_tilde
                    * (
                        D
                        * (1 + gamma * pole_disp * tc)
                        * (
                            dtcdr * Nr2
                            + (dNr2dr - dtcdr * gamma * Nr2 * pole_disp) * tc
                            + dNr2dr * gamma * pole_disp * tc ^ 2
                        )
                        - Nr2
                        * tc
                        * (
                            dtcdr * gamma ^ 2 * Nr2 * pole_disp ^ 2 * tc
                            + dDdr * (1 + gamma * pole_disp * tc) ^ 2
                        )
                    )
                )
            )
        )
        / (
            (
                1
                + 2 * gamma * pole_disp * tc
                + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
            )
            * (
                gamma ^ 2 * Nr2 ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                + 2
                * D
                * gamma
                * Nr2
                * tc
                * (pole_disp + gamma * (omega_tilde ^ 2 + pole_disp ^ 2) * tc)
                + D ^ 2
                * (
                    1
                    + 2 * gamma * pole_disp * tc
                    + gamma ^ 2 * (omega_tilde ^ 2 + pole_disp ^ 2) * tc ^ 2
                )
            )
        )
    ) / ((omega_tilde ^ 2 + pole_disp ^ 2) ^ 2 * r ^ 2)
    return C1R, C1I, C0R, C0I, PsiR, PsiI
end
