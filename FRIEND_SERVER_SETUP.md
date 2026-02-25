# Friend Setup: Run Codex + Analysis on Institute Server (250GB data)

This is the exact walkthrough to run the analysis on a new computer.

## 1) Local prerequisites (friend laptop)
Install:
- `ssh`, `rsync`, `git`, `tmux`
- VPN client used by your institute

macOS (Homebrew):
```bash
brew install git rsync tmux
```

## 2) SSH key setup (if needed)
```bash
ls ~/.ssh/id_ed25519.pub || ssh-keygen -t ed25519 -C "friend@laptop"
```
Add `~/.ssh/id_ed25519.pub` to institute account authorized keys.

## 3) SSH host alias
Create/edit `~/.ssh/config` and add:
```sshconfig
Host institute-bio
  HostName <SERVER_FQDN_OR_IP>
  User <USERNAME>
  IdentityFile ~/.ssh/id_ed25519
  ServerAliveInterval 30
  ServerAliveCountMax 120
```

## 4) Connect VPN, then verify SSH
```bash
ssh institute-bio 'hostname; whoami'
```

## 5) Copy project to server
From local machine:
```bash
rsync -avz /path/to/local/bio/ institute-bio:~/projects/crispr_analysis/
```

## 6) Start persistent server shell
```bash
ssh institute-bio
tmux new -s crispr
cd ~/projects/crispr_analysis
```

Detach safely with `Ctrl-b d`, resume with:
```bash
tmux attach -t crispr
```

## 7) Set full dataset path
Edit server file:
`~/projects/crispr_analysis/results1/sceptre_realdata/config.yaml`

Update:
- `data.h5ad_path` -> full 250GB dataset path on server
- `data.results_root` -> desired output directory on server storage

## 8) Ensure Hallmark GMT exists on server
```bash
mkdir -p ~/projects/crispr_analysis/results1/sceptre_realdata/pathways
curl -L --fail https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt \
  -o ~/projects/crispr_analysis/results1/sceptre_realdata/pathways/h.all.v2024.1.Hs.symbols.gmt
```

## 9) Run the R-canonical full pipeline
```bash
cd ~/projects/crispr_analysis/results1/sceptre_realdata
bash run_all_r_canonical.sh
```

## 10) Monitor progress
```bash
ls -lt ~/projects/crispr_analysis/results1/sceptre_realdata/sceptre | head
ls -lt ~/projects/crispr_analysis/results1/sceptre_realdata/scanpy | head
ls -lt ~/projects/crispr_analysis/results1/sceptre_realdata/figures | head
```

## 11) Pull results back to laptop
```bash
rsync -avz institute-bio:~/projects/crispr_analysis/results1/sceptre_realdata/ ./sceptre_realdata_results/
```

## Important
- Do not run heavy analysis through Samba mounts.
- Samba is OK only for browsing/copying final outputs.
- Optional: copy `ssh_config.institute-bio.template` into `~/.ssh/config` and fill placeholders.
