#!/bin/bash

eval "$(ssh-agent -s)"

ssh-add ~/.ssh/id_rsa

git fetch origin
git reset --hard origin/main