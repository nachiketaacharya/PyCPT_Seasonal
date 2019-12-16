      Module cpt_lapack
!
      Interface
         Subroutine dgesdd ( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT,   &
     &                      WORK, LWORK, IWORK, INFO )
            Character          JOBZ
            Integer            INFO, LDA, LDU, LDVT, LWORK, M, N
            Integer            IWORK( * )
            Double Precision   A( LDA, * ), S( * ), U( LDU, * ),        &
     &                         VT( LDVT, * ), WORK( * )
         End Subroutine dgesdd
!
         Subroutine sgesdd ( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT,   &
     &                      WORK, LWORK, IWORK, INFO )
            Character          JOBZ
            Integer            INFO, LDA, LDU, LDVT, LWORK, M, N
            Integer            IWORK( * )
            Real               A( LDA, * ), S( * ), U( LDU, * ),        &
     &                         VT( LDVT, * ), WORK( * )
         End Subroutine sgesdd
      End Interface
!
      End Module cpt_lapack
