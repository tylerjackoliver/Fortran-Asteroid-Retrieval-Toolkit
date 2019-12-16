        module precision_kinds

            integer, parameter :: &

                sp = kind(1.0),                              &
                dp = selected_real_kind(2*precision(1.0_sp)),&
                qp = selected_real_kind(4*precision(1.0_sp))

        end module precision_kinds
