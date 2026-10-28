"""Check repository documentation links/anchors and optionally run tutorials.

Uses only the standard library. External URLs are intentionally not requested
in CI; bibliographic metadata is reviewed against primary sources separately.
"""
import argparse
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path
import re
import subprocess
import sys
from urllib.parse import unquote, urlsplit


class Tags(HTMLParser):
    def __init__(self):
        super().__init__(); self.links=[]; self.ids=set()

    def handle_starttag(self, tag, attrs):
        data=dict(attrs)
        if data.get('id'):
            self.ids.add(data['id'])
        for key in ('href','src'):
            if data.get(key):
                self.links.append(data[key])


def parse(path):
    text=path.read_text(encoding='utf-8')
    # Code snippets contain code, not rendered Markdown headings/links.
    body=re.sub(r'^```[^\n]*\n.*?^```[^\n]*$', '', text, flags=re.M|re.S)
    tags=Tags(); tags.feed(body)
    ids=set(tags.ids); counts=Counter()
    for heading in re.findall(r'^#{1,6}\s+(.+?)\s*#*$',body,flags=re.M):
        heading=re.sub(r'\[([^]]+)\]\([^)]*\)',r'\1',heading)
        slug=re.sub(r'[^\w\- ]','',heading.lower()).replace(' ','-')
        index=counts[slug];counts[slug]+=1
        ids.add(slug+(f'-{index}' if index else ''))
    links=tags.links+re.findall(r'!?\[[^\]\n]*\]\(([^)\n]+)\)',body)
    links+=re.findall(r'^\[[^]]+\]:\s*(\S+)',body,flags=re.M)
    return text, ids, links


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-examples',action='store_true')
    args=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    paths=sorted(set(root.glob('*.md'))|set((root/'docs').rglob('*.md'))|
                 set((root/'examples').rglob('*.md')))
    parsed={path:parse(path) for path in paths}
    errors=[]
    for source,(_,_,links) in list(parsed.items()):
        for raw in links:
            # Optional Markdown title follows the URL; spaces in a local
            # filename should use the supported <...> form.
            url=raw[1:raw.index('>')] if raw.startswith('<') and '>' in raw else raw.split()[0]
            parts=urlsplit(url)
            if parts.scheme or parts.netloc:
                continue
            target=(source.parent/unquote(parts.path)).resolve() if parts.path else source
            if not target.exists():
                errors.append(f'{source.relative_to(root)}: missing {url}')
            elif parts.fragment and target.suffix=='.md':
                if target not in parsed:
                    parsed[target]=parse(target)
                if unquote(parts.fragment) not in parsed[target][1]:
                    errors.append(f'{source.relative_to(root)}: missing anchor {url}')
    keys=re.findall(r'@\w+\{([^,]+),',(root/'references.bib').read_text(encoding='utf-8'))
    if len(keys)!=len(set(keys)):
        errors.append('references.bib: duplicate keys')
    ref_ids=parsed[root/'docs/references.md'][1]
    errors.extend(f'references.bib: no annotated entry for {key}' for key in keys if key not in ref_ids)
    scripts=sorted((root/'examples/howto').glob('[0-9][0-9]_*.py'))
    howto=parsed[root/'HOWTO.md'][0]
    errors.extend(f'HOWTO.md: tutorial not linked: {script.name}' for script in scripts
                  if script.relative_to(root).as_posix() not in howto)
    if errors:
        raise SystemExit('\n'.join(errors))
    print(f'Checked {len(paths)} Markdown files, {len(keys)} bibliography entries and {len(scripts)} linked tutorials.',flush=True)
    if args.run_examples:
        for script in scripts:
            choices=['jpa','jtwpa','rf-squid','snail'] if script.name.startswith('08_') else [None]
            for family in choices:
                command=[sys.executable,str(script)]
                if family:
                    command+=['--family',family]
                print('Running',script.relative_to(root),family or '',flush=True)
                subprocess.run(command,cwd=root,check=True,timeout=180)
        print('All HOWTO examples completed; JoSIM was dry-run only.')


if __name__=='__main__':
    main()
