import sys
import click
import logging
from authormaps.authorinfo import AuthorData
from authormaps.author_network import AuthorNetwork
from authormaps.networkutils import *


@click.group()
def cli():
    pass


@cli.command(name="get_publications")
@click.argument("author_name", type=str)
def get_publications(author_name):
    author = AuthorData(author_name)
    publications = author.get_list_of_publications()
    if publications:
        print(publications)
    else:
        print(f"No publications for '{author_name}'")


@cli.command(name="get_coauthors")
@click.argument("author_name", type=str)
def get_coauthors(author_name):
    author = AuthorData(author_name)
    coauthors = author.get_list_of_coauthors_from_list_of_publications()
    if coauthors:
        print(coauthors)
    else:
        print(f"No coauthors for '{author_name}'")


@cli.command(name="get_common_publications")
@click.argument("author_name_1", type=str)
@click.argument("author_name_2", type=str)
def get_common_publications(author_name_1, author_name_2):
    publications = get_publications_for_pair((author_name_1, author_name_2))
    if publications:
        print(publications)
    else:
        print(f"No common publications for '{author_name_1}' and '{author_name_2}'")


@cli.command(name="get_common_publications_for_all_coauthors")
@click.argument("author_name", type=str)
def get_common_publications_for_all_coauthors(author_name):
    author = AuthorData(author_name)
    coauthors = author.get_list_of_coauthors_from_list_of_publications()
    if coauthors:
        shared_publications = get_shared_publications_for_authors_pairs(coauthors)
        if shared_publications:
            print(shared_publications)
        else:
            print(f"No common publications for '{author_name}'s coauthors")
    else:
        print(f"No coathors for '{author_name}'")


@cli.command(name="get_common_publications_counts_for_all_coauthors")
@click.argument("author_name", type=str)
def get_common_publications_counts_for_all_coauthors(author_name):
    author = AuthorData(author_name)
    coauthors = author.get_list_of_coauthors_from_list_of_publications()
    if coauthors:
        shared_publications_counts = count_shared_publications(coauthors)
        if shared_publications_counts:
            print(shared_publications_counts)
        else:
            print(f"No common publications for '{author_name}'s coauthors")
    else:
        print(f"No coathors for '{author_name}'")


@cli.command(name="generate_graph")
@click.argument("author_name", type=str)
@click.argument("output_file", type=click.Path())
@click.option("-f", "--format", type=click.Choice(["png", "jpg", "svg", "pdf"]), default="png")
def generate_graph(author_name, output_file, format):
    author = AuthorData(author_name)
    coauthors = author.get_list_of_coauthors_from_list_of_publications()
    if coauthors:
        shared_publications_counts = count_shared_publications(coauthors)
        if shared_publications_counts:
            network = AuthorNetwork(shared_publications_counts)
            network.save_graph(output_format=format, filename=output_file)
        else:
            print(f"No common publications for '{author_name}'s coauthors")
    else:
        print(f"No coathors for '{author_name}'")


def main():
    cli()


if __name__ == '__main__':
    main()
