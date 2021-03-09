
execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    TIMEOUT 1
    OUTPUT_VARIABLE BRANCH_NAME
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    )

execute_process(
    COMMAND git rev-list HEAD -n 1
    TIMEOUT 1
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    )

string(TIMESTAMP OURTIME %Y-%m-%d\ %H:%M:%S)
string(TIMESTAMP UTCTIME %Y-%m-%d\ %H:%M:%S UTC)
