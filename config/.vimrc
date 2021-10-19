set expandtab
set tabstop=4
set smartindent
set shiftwidth=4
set complete-=i
set number
syntax on
filetype plugin indent on
" Plugins will be downloaded under the specified directory.
call plug#begin('~/.vim/bundle')
" Declare the list of plugins.
Plug 'tpope/vim-sensible'
Plug 'junegunn/seoul256.vim'
Plug 'iamcco/mathjax-support-for-mkdp'
Plug 'iamcco/markdown-preview.vim'
Plug 'pechorin/any-jump.vim'
" List ends here. Plugins become visible to Vim after this call.
call plug#end()
autocmd BufNewFile,BufRead *.md set filetype=markdown
