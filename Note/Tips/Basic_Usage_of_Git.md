## Motivation
1. How to submit code in gitlab/github.
2. How to record the modification.
3. How to track submission.
4. How to synchronize the code.

## Scenario Demo

### Step 1: Grab source code from public repo to local
1. Prepare 2 local working-dirs, including:
   * Production Dir
     * Purpose: run production data
   * Development Dir
     * Purpose: for development
   * Git Command line
   ```
   e.g.
   gitRepo="git@10.133.130.114:lxwgcool/RNA_Auto_Launcher.git"
   $git clone ${gitRepo}
   ```

2. For Production Dir
   * Use the code to run prodicton flowcells
   * For all standard analysis, which need to shared by others
   * Never edit any files manually in this folder   

3. For Development Dir
   * Use this code for development, including
     * Debug
     * New Feature
     * Code investigation 
     * Any other development works

### Step 2: Create a gitlab/github issue (ticket)
1. Descrbe the detail of the requirement and record it
2. remember the number of this issue. 

### Step 3: Development
1. Impelemtn the development work, such as debug, add new feature and etc. 

### Step 4: 






