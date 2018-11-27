import re

def edit_qiime_params(param):

    if not re.match("\-\-\w\-", param):
        return False

    return re.sub("\-\-\w\-", "", param)

def edit_qiime_types(qtype):

    if not re.match(pattern="\w+(?:\[\w+\])?",string=qtype):
        raise Exception("type {type} is not a valid format".format(type=qtype))

    return re.sub(pattern="\[",
                  repl=".",
                  string=re.sub(pattern="\]",
                                repl="",
                                string=qtype))