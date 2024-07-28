Here I will be describing the method to the vault.

# General

# Specifics

## Packages

For packages, there is certain information I always include. Some are:
1. **Metadata**: In the frontmatter, I indicate that this page describes a package and I note if the packages is being maintained or not.
2. **Installation**: If you get a new device or the HPC moves your home directory and messes up your conda environments (happened to me). 
3. **Documentation Links**: Packages get updated, it always good to revisit their documentation. For example, the developer might add a feature that you actually really need or didn’t know you need.
4. **Usage examples**: I like to write done commands that I use but can’t seem to remember. It’s also a good place to copy and paste from. 
5. **Errors**: Some times a random error shows up and you don’t see it again for months. When the error shows up again, don’t waste your time googling how to fix it again. 

## Research papers

A really useful part of this vault is how to organize and use research papers within the vault. This was implemented after coming across [Christian B. B. Houmann’s tutorial](https://bagerbach.com/blog/how-i-read-research-papers-with-obsidian-and-zotero) (an important note is that I implemented the workflow before it was revamped in 2023). This requires the citation manager [Zotero](https://www.zotero.org) in conjunction with [Better BibTeX (BBT)](https://retorque.re/zotero-better-bibtex/). 

First, [install BBT](https://retorque.re/zotero-better-bibtex/installation/) then export your Zotero library.

After you set up Zotero, you will need to point to your exported Zotero library under the citation plugin settings:
![[citation plugin settings.png|500]]

Once you have done that, to add a paper to your vault, press the hotkey command + P and type in “Add paper” in the pop up window:

![[QuickAdd add paper.png|600]]

After selecting “QuickAdd: Add Paper”, start searching for the paper you want to add by typing a word from the paper’s title:

![[QuickAdd add paper search.png|600]]
A new page with the paper’s citation key will be created in your vault under the folder Research/Papers.