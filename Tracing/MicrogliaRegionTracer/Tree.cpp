#include "Tree.h"

Tree::Tree()
{
	this->root = NULL;
}

Tree::~Tree()
{
	//TODO: Delete all the nodes in the tree, otherwise big memory leak!
}

void Tree::SetRoot(Node *root)
{
	member_nodes.push_back(root);
	this->root = root;
}

void Tree::AddNode(Node* node, Node* parent)
{
	//The node itself is responsible for keeping track of its parent, as well as the parent keeping track of the node
	
	if (parent == NULL)
		throw std::runtime_error("AddNode needs a parent");

	member_nodes.push_back(node);
}

std::vector<Node*> Tree::GetMemberNodes()
{
	return member_nodes;
}

Node* Tree::getRoot()
{
	return root;
}

bool Tree::RemoveNode(Node* node)
{
	std::vector< Node* >::iterator member_nodes_iter;
	for (member_nodes_iter = member_nodes.begin(); member_nodes_iter != member_nodes.end(); ++member_nodes_iter)
	{
		Node* member_node = *member_nodes_iter;

		if (member_node == node)
		{
			member_nodes.erase(member_nodes_iter);
			return true;
		}
	}

	return false;
}