{% extends 'markdown.tpl' %}

{% block input %}
<pre><code class="language-python">{{ cell.source }}</code>
</pre>
{% endblock input %}

{% block stream %}
{% if output.name == "stdout" %}
<pre><code class="language-">{{ output.text.strip() | escape }}</code>
</pre>
{% endif %}
{% endblock stream %}

{% block data_text scoped %}
<pre><code class="language-">{{ output.data['text/plain'].strip() | escape }}</code>
</pre>
{% endblock data_text %}

