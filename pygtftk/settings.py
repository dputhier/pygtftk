# ---------------------------------------------------------------
# The script for bash completion
# ---------------------------------------------------------------


COMPLETION_SCRIPT = """
# To activate completion of sub-commands
# Add this script to your ~/.bashrc
 _gtftk()
 {{
     local cur;
     COMPREPLY=();
     cur=${{COMP_WORDS[COMP_CWORD]}};
     if [ ${{COMP_CWORD}} == 1 ]; then
     COMPREPLY=( $( compgen -W '{p}' -- ${{cur}}) )
     else
     COMPREPLY=( $( compgen -f -o plusdirs -- ${{cur}} ) )
     fi
 }}
 complete -o nospace -F _gtftk pygtftk
"""


def get_completion_script():
    import pygtftk.cmd_manager
    prg_list = pygtftk.cmd_manager.CmdManager.cmd_obj_list.keys()
    prg_str = " ".join(sorted(prg_list))
    return COMPLETION_SCRIPT.format(p=prg_str)
