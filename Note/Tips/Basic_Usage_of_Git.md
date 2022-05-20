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
1. Synchronize you local repo with public repo first.
```
git pull
```
3. Implement the development work, such as debug, adding new features and etc. 
4. Remember keep recording some important messages (if there are) in the gitlab/github issue (ticket)
5. All modifications should base on the source code in "Development Dir"
.
### Step 4: Synchronize the new modification with your local git repo 
1. Enter into your "Development Dir"
2. Synchronize your local git repo with public repo
```
git pull
```
3. Synchronize your local code with local git repo
```
git add -A
```
4. make comments 
    * Notice
       * Format, including caption, main body and ref
       * do not write too long for each line
       * make sure you linked this commit to the correct ticket.
       * Save it as Dos format
    * Command line and example 
```
Command line:
git commit

Comments Example: 
commit f75c4ad7331f525932e0c4cd0689a221f4b07068
Author: lxwgcool <lxwgcool@gmail.com>
Date:   Fri May 20 14:43:36 2022 -0400

    Demo: update multi-files
    Details:
    1: Try to modify 3 differnet files including
       1) Deploy_Auto_Launcher_in_HPC.MD
       2) SourceCode/AutoMinRNAPipelineLauncher.py
       3) SourceCode/global_config_bash.rc
    2: Simply comments, not impact any real functions of the code.
    
    ref #11
```

### Step 5: Synchronize your local git repo with public git repo
1. Command line 
```
git push
```

### Step 6: Synchronize the code in your production dir
1. Go to your production dir and synchronize the code with public repo
```
git pull
```

### Step 7: Use the source code in production dir to run real data
Just do it and enjoy!
