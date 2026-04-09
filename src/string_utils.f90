MODULE string_utils
  IMPLICIT NONE
CONTAINS

  ! Converts an integer to a trimmed string
  FUNCTION itoa(i) RESULT(str)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    CHARACTER(LEN=20) :: str
    WRITE(str, '(I0)') i  ! I0: no leading spaces
  END FUNCTION itoa

  ! Replaces all occurrences of one character with another in a string
  SUBROUTINE replace_char(str, old_char, new_char)
    IMPLICIT NONE
    CHARACTER(*), INTENT(INOUT) :: str
    CHARACTER(1), INTENT(IN)    :: old_char, new_char
    INTEGER :: i

    DO i = 1, LEN_TRIM(str)
      IF (str(i:i) == old_char) str(i:i) = new_char
    END DO
  END SUBROUTINE replace_char

END MODULE string_utils
