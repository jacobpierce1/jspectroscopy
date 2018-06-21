git filter-branch -f --env-filter \
"GIT_AUTHOR_NAME='Jacob'; GIT_AUTHOR_EMAIL='jacobpierce@uchicago.edu'; \
GIT_COMMITTER_NAME='Jacob'; GIT_COMMITTER_EMAIL='jacobpierce@uchicago.edu';" HEAD