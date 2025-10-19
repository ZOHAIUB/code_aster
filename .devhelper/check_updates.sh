#!/bin/bash

set_prefix() {
    local this=$(readlink -n -f "$1")
    devhelper=$(dirname "${this}")
    root=$(dirname $(dirname "${devhelper}"))
}

set_prefix "${0}"

test -z "${CI_PROJECT_ID}" || exit 0
test -z "${CI_CODEASTER}" || exit 0

mark="${devhelper}/.last-check-updates"
markbase=$(basename "${mark}")

NB_DAYS=5

update_mark()
{
    echo "updating mark..."
    date > "${mark}"
}

check_updates_main()
{
    local found
    printf "checking for tools updates... "
    if [ ! -f "${mark}" ]; then
        echo "never run, updating..."
        do_updates
    else
        found=$( find "${devhelper}" -mtime +${NB_DAYS} -name "${markbase}" )
        if [ ! -z "${found}" ] || [ ! -z "${FORCE_CHECK_UPDATES}" ]; then
            echo "not checked in last 10 days, updating..."
            do_updates
        else
            echo "not needed"
        fi
    fi
    if [ ! -e .git/hooks/pre-commit ]; then
        echo "Git hooks not installed, installing..."
        do_install_env
    fi
}

do_updates()
{
    update_repo src
    update_repo devtools
    update_repo data
    update_repo validation
    do_install_env
    update_mark
}

update_repo()
{
    local repo="$1"
    cd "${root}"
    if [ ! -d ${repo} ]; then
        return
    fi
    printf "updating ${repo}... "
    current_branch=$(git -C ${repo} rev-parse --abbrev-ref HEAD)
    if [ "${repo}" != "src" ] && [ "${current_branch}" = "main" ]; then
        printf "updating 'main'... "
        command=pull
    else
        printf "fetching 'main'... "
        command=fetch
    fi
    git -C ${repo} ${command} origin main > /dev/null 2>&1
    echo "completed"
    if [ "${repo}" != "src" ] && [ "${current_branch}" != "main" ]; then
        echo "WARNING: ${repo} is not on 'main'"
    fi
}

do_install_env()
{
    if [ ! -d "${root}/devtools" ]; then
        echo "ERROR: devtools not found!"
        return
    fi
    cd "${root}/src"
    export ASTER_CHECK_UPDATES=1
    "${root}/devtools/bin/install_env"
}

check_updates_main "${@}"
