all: post_root post_root_nolen

post_root: 
	gcc post_root.c -o post_root

post_root_nolen: 
	gcc post_root_nolen.c -o post_root_nolen

clean:
	rm post_root post_root_nolen
